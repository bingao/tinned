#include <iostream>

#include <cstddef>
#include <string>

#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/constants.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/matrices/identity_matrix.h>
#include <symengine/matrices/diagonal_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_exception.h>

#include "Tinned/TwoLevelAtom.hpp"

namespace Tinned
{
    TwoLevelOperator::TwoLevelOperator(
        const std::pair<SymEngine::RCP<const OneElecOperator>,
                        SymEngine::RCP<const SymEngine::MatrixExpr>>& H0,
        const std::map<SymEngine::RCP<const OneElecOperator>,
                       SymEngine::RCP<const SymEngine::MatrixExpr>,
                       SymEngine::RCPBasicKeyLess>& V,
        const std::pair<SymEngine::RCP<const OneElecDensity>,
                        SymEngine::RCP<const SymEngine::MatrixExpr>>& rho0
    ) : max_order_(7), H0_(H0), V_(V)
    {
        // A density matrix is Hermitian and has a trace one
        if (SymEngine::neq(
            *rho0.second, *SymEngine::conjugate_matrix(SymEngine::transpose(rho0.second))
        )) throw SymEngine::SymEngineException(
                "Density matrix must be Hermitian: " + stringify(rho0.second)
            );
        if (!is_zero_quantity(SymEngine::add(
            SymEngine::minus_one, SymEngine::trace(rho0.second)
        ))) throw SymEngine::SymEngineException(
                "Density matrix must have a trace one: " + stringify(rho0.second)
            );
        rho0_ = rho0;
        auto op = SymEngine::rcp_dynamic_cast<const SymEngine::DiagonalMatrix>(H0.second);
        // Here, we do not check if all matrices should have the same dimension
        dim_operator_ = op->get_container().size();
        for (std::size_t i=0; i<dim_operator_; ++i) {
            for (std::size_t j=0; j<dim_operator_; ++j) {
                if (i==j) {
                    omega_.push_back(SymEngine::zero);
                }
                else {
                    omega_.push_back(SymEngine::sub(op->get(i), op->get(j)));
                }
            }
        }
        SymEngine::vec_basic values;
        for (const auto& oper: V_) {
            auto dependencies = oper.first->get_dependencies();
            if (dependencies.size()!=1) throw SymEngine::SymEngineException(
                "Each field operator should depend only on one perturbation: "
                + stringify(oper.first)
            );
            if (perturbations_.find(dependencies.begin()->first)!=perturbations_.end())
                throw SymEngine::SymEngineException(
                    "Field operators should depend on different perturbations: "
                    + stringify(oper.first)
                );
            perturbations_.insert(dependencies.begin()->first);
            //Elements of V's can be symbols that we cannot check the hermicity
            //// All V's must be Hermitian
            //if (SymEngine::neq(
            //    *oper.second,
            //    *SymEngine::conjugate_matrix(SymEngine::transpose(oper.second))
            //)) throw SymEngine::SymEngineException(
            //        "Each field operator must be Hermitian: " + stringify(oper.second)
            //    );
            values.push_back(oper.second);
        }
        // Each pair of V's must commute
        for (std::size_t i=1; i<values.size(); ++i) {
            for (std::size_t j=0; j<i; ++j) {
                auto val_ij = SymEngine::matrix_mul({values[i], values[j]});
                auto val_ji = SymEngine::matrix_mul({values[j], values[i]});
                if (SymEngine::neq(*val_ij, *val_ji)) throw SymEngine::SymEngineException(
                    "Field operators: " + stringify(values[i]) + " and "
                    + stringify(values[j]) + " do not commute"
                );
            }
        }
        // Store unperturbed density matrix
        rho_all_derivatives_[0] = DensityDerivative({
            std::make_pair(SymEngine::multiset_basic({}), rho0.second)
        });
    }

    SymEngine::vec_basic TwoLevelOperator::get_values(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& A
    ) const
    {
        if (SymEngine::is_a_sub<const SymEngine::ZeroMatrix>(*A)) {
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ZeroMatrix>(A);
            if (SymEngine::neq(*op->nrows(), *SymEngine::integer(dim_operator_)) ||
                SymEngine::neq(*op->ncols(), *SymEngine::integer(dim_operator_)))
                throw SymEngine::SymEngineException(
                    "Incorrect dimension (" + std::to_string(dim_operator_)
                    + ") of matrix: " + stringify(A)
                );
            SymEngine::vec_basic values;
            for (std::size_t i=0; i<dim_operator_*dim_operator_; ++i)
                values.push_back(SymEngine::zero);
            return values;
        }
        else if (SymEngine::is_a_sub<const SymEngine::IdentityMatrix>(*A)) {
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::IdentityMatrix>(A);
            if (SymEngine::neq(*op->size(), *SymEngine::integer(dim_operator_)))
                throw SymEngine::SymEngineException(
                    "Incorrect dimension (" + std::to_string(dim_operator_)
                    + ") of matrix: " + stringify(A)
                );
            SymEngine::vec_basic values;
            for (std::size_t i=0; i<dim_operator_; ++i) {
                for (std::size_t j=0; j<dim_operator_; ++j) {
                    if (i==j) {
                        values.push_back(SymEngine::one);
                    }
                    else {
                        values.push_back(SymEngine::zero);
                    }
                }
            }
            return values;
        }
        else if (SymEngine::is_a_sub<const SymEngine::DiagonalMatrix>(*A)) {
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::DiagonalMatrix>(A);
            auto container = op->get_container();
            if (container.size()!=dim_operator_)
                throw SymEngine::SymEngineException(
                    "Incorrect dimension (" + std::to_string(dim_operator_)
                    + ") of matrix: " + stringify(A)
                );
            SymEngine::vec_basic values;
            for (std::size_t i=0; i<container.size(); ++i) {
                for (std::size_t j=0; j<container.size(); ++j) {
                    if (i==j) {
                        values.push_back(container[i]);
                    }
                    else {
                        values.push_back(SymEngine::zero);
                    }
                }
            }
            return values;
        }
        else if (SymEngine::is_a_sub<const SymEngine::ImmutableDenseMatrix>(*A)) {
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
            if (op->nrows()!=dim_operator_ || op->ncols()!=dim_operator_)
                throw SymEngine::SymEngineException(
                    "Incorrect dimension (" + std::to_string(dim_operator_)
                    + ") of matrix: " + stringify(A)
                );
            return op->get_values();
        }
        else {
            throw SymEngine::SymEngineException("Unsupported matrix: "+stringify(A));
        }
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_hermitian_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A)
    {
        auto conj = SymEngine::conjugate_matrix(A);
        return SymEngine::transpose(conj);
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_1el_density(const OneElecDensity& x)
    {
        auto derivatives = x.get_derivatives();
        if (derivatives.size()>max_order_) SymEngine::SymEngineException(
            "Invalid order: " + std::to_string(derivatives.size())
            + "(" + std::to_string(max_order_) + ")"
        );
        switch (derivatives.size()) {
            case 0:
                return rho0_.second;
            default:
                if (derivatives.size()>rho_all_derivatives_.size()-1) {
                    for (unsigned int order=rho_all_derivatives_.size()-1;
                         order<derivatives.size();
                         ++order) {
                        DensityDerivative rho_derivatives;
                        auto pert_permutation = PertPermutation(order+1, perturbations_);
                        bool remaining = true;
                        do {
                            auto permut_derivatives
                                = pert_permutation.get_derivatives(remaining);
                            auto freq_sum = get_frequency_sum(permut_derivatives);
                            // The last perturbation is the one for the
                            // external field, and we need to find the value of
                            // its corresponding operator
                            SymEngine::RCP<const SymEngine::MatrixExpr> oper_field;
                            bool not_found = true;
                            for (const auto& oper: V_) {
                                if (SymEngine::eq(
                                    *permut_derivatives.back(),
                                    *oper.first->get_dependencies().begin()->first
                                )) {
                                    oper_field = oper.second;
                                    not_found = false;
                                    break;
                                }
                            }
                            if (not_found) throw SymEngine::SymEngineException(
                                "Invalid perturbation for the external field "
                                + stringify(permut_derivatives.back())
                            );
                            // Check if the higher order derivatives of the
                            // density matrix already computed partially
                            auto higher_derivatives = SymEngine::multiset_basic(
                                permut_derivatives.begin(), permut_derivatives.end()
                            );
                            int idx_higher = -1;
                            for (std::size_t i=0; i<rho_derivatives.size(); ++i)
                                if (SymEngine::unified_eq(rho_derivatives[i].first, higher_derivatives)) {
                                    idx_higher = i;
                                    break;
                                }
                            // Try to find the matched lower order derivatives
                            // of the density matrix
                            not_found = true;
                            auto lower_derivatives = SymEngine::multiset_basic(
                                permut_derivatives.begin(), permut_derivatives.end()-1
                            );
                            for (const auto& rho_lower: rho_all_derivatives_[order]) {
                                if (SymEngine::unified_eq(rho_lower.first, lower_derivatives)) {
                                    auto val_V_rho = get_values(SymEngine::matrix_mul({
                                        oper_field, rho_lower.second
                                    }));
                                    auto val_rho_V = get_values(SymEngine::matrix_mul({
                                        rho_lower.second, oper_field
                                    }));
                                    SymEngine::vec_basic val_rho;
                                    // Compute higher order derivatives of
                                    // the density matrix
                                    for (std::size_t i=0; i<omega_.size(); ++i) {
                                        val_rho.push_back(
                                            SymEngine::div(
                                                SymEngine::sub(val_V_rho[i], val_rho_V[i]),
                                                SymEngine::sub(freq_sum, omega_[i])
                                            )
                                        );
                                    }
                                    if (idx_higher>=0) {
                                        rho_derivatives[idx_higher].second
                                            = SymEngine::matrix_add({
                                                  rho_derivatives[idx_higher].second,
                                                  SymEngine::immutable_dense_matrix(
                                                      dim_operator_, dim_operator_, val_rho
                                                  )
                                              });
                                    }
                                    else {
                                        rho_derivatives.push_back(
                                            std::make_pair(
                                                higher_derivatives,
                                                SymEngine::immutable_dense_matrix(
                                                    dim_operator_, dim_operator_, val_rho
                                                )
                                            )
                                        );
                                    }
                                    not_found = false;
                                    break;
                                }
                            }
                            if (not_found) {
                                std::string str_lower_derivatives;
                                for (const auto& p: permut_derivatives)
                                   str_lower_derivatives += stringify(p) + ",";
                                throw SymEngine::SymEngineException(
                                    "Invalid lower order derivatives: "
                                    + str_lower_derivatives
                                    + " of the density matrix: "
                                    + stringify(x)
                                );
                            }
                        } while (remaining);
                        rho_all_derivatives_[order+1] = rho_derivatives;
                    }
                }
                for (const auto& rho_derivatives: rho_all_derivatives_[derivatives.size()])
                    if (SymEngine::unified_eq(derivatives, rho_derivatives.first))
                        return rho_derivatives.second;
                throw SymEngine::SymEngineException(
                    "Density matrix with invalid derivatives: " + stringify(x)
                );
        }
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_1el_operator(const OneElecOperator& x)
    {
        if (x.get_name()==H0_.first->get_name()) {
            auto derivatives = x.get_derivatives();
            if (derivatives.empty()) {
                return H0_.second;
            }
            else {
                return SymEngine::zero_matrix(
                    SymEngine::integer(dim_operator_), SymEngine::integer(dim_operator_)
                );
            }
        }
        else {
            for (const auto& oper: V_) {
                if (x.get_name()==oper.first->get_name()) {
                    auto derivatives = x.get_derivatives();
                    if (derivatives.empty()) {
                        return SymEngine::matrix_mul({
                            oper.first->get_dependencies().begin()->first,
                            oper.second
                        });
                    }
                    else if (derivatives.size()==1) {
                        return oper.second;
                    }
                    else {
                        return SymEngine::zero_matrix(
                            SymEngine::integer(dim_operator_),
                            SymEngine::integer(dim_operator_)
                        );
                    }
                }
            }
        }
        throw SymEngine::SymEngineException("Invalid operator " + stringify(x));
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_temporum_operator(const TemporumOperator& x)
    {
        result_ = apply(x.get_target());
        eval_oper_scale(x.get_frequency(), result_);
        return result_;
    }

    SymEngine::RCP<const SymEngine::MatrixExpr> TwoLevelOperator::eval_conjugate_matrix(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& A
    )
    {
        return SymEngine::conjugate_matrix(A);
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A)
    {
        return SymEngine::transpose(A);
    }

    void TwoLevelOperator::eval_oper_addition(
        SymEngine::RCP<const SymEngine::MatrixExpr>& A,
        const SymEngine::RCP<const SymEngine::MatrixExpr>& B
    )
    {
        auto val_A = get_values(A);
        auto val_B = get_values(B);
        SymEngine::vec_basic values;
        for (std::size_t i=0; i<val_A.size(); ++i)
            values.push_back(SymEngine::add(val_A[i], val_B[i]));
        A = SymEngine::immutable_dense_matrix(dim_operator_, dim_operator_, values);
    }

    SymEngine::RCP<const SymEngine::MatrixExpr> TwoLevelOperator::eval_oper_multiplication(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& A,
        const SymEngine::RCP<const SymEngine::MatrixExpr>& B
    )
    {
        return SymEngine::matrix_mul({A, B});
    }

    void TwoLevelOperator::eval_oper_scale(
        const SymEngine::RCP<const SymEngine::Basic>& scalar,
        SymEngine::RCP<const SymEngine::MatrixExpr>& A
    )
    {
        auto val_A = get_values(A);
        SymEngine::vec_basic values;
        for (std::size_t i=0; i<val_A.size(); ++i)
            values.push_back(SymEngine::mul(scalar, val_A[i]));
        A = SymEngine::immutable_dense_matrix(dim_operator_, dim_operator_, values);
    }

    TwoLevelFunction::TwoLevelFunction(
        const std::pair<SymEngine::RCP<const OneElecOperator>,
                        SymEngine::RCP<const SymEngine::MatrixExpr>>& H0,
        const std::map<SymEngine::RCP<const OneElecOperator>,
                       SymEngine::RCP<const SymEngine::MatrixExpr>,
                       SymEngine::RCPBasicKeyLess>& V,
        const std::pair<SymEngine::RCP<const OneElecDensity>,
                        SymEngine::RCP<const SymEngine::MatrixExpr>>& rho0
    ) : FunctionEvaluator<SymEngine::RCP<const SymEngine::Basic>,
                          SymEngine::RCP<const SymEngine::MatrixExpr>>(
            std::make_shared<TwoLevelOperator>(H0, V, rho0)
        ) {}

    SymEngine::RCP<const SymEngine::Basic> TwoLevelFunction::eval_trace(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& A
    )
    {
        return SymEngine::trace(A);
    }

    void TwoLevelFunction::eval_fun_addition(
        SymEngine::RCP<const SymEngine::Basic>& f,
        const SymEngine::RCP<const SymEngine::Basic>& g
    )
    {
        f = SymEngine::add(f, g);
    }

    void TwoLevelFunction::eval_fun_scale(
        const SymEngine::RCP<const SymEngine::Basic>& scalar,
        SymEngine::RCP<const SymEngine::Basic>& f
    )
    {
        f = SymEngine::mul(scalar, f);
    }
}
