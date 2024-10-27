#include <cstddef>
#include <iterator>
#include <string>

#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/integer.h>
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
#include <symengine/simplify.h>

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
    ) : max_order_(12), H0_(H0), V_(V)
    {
        // Compute transition angular frequencies
        auto op = SymEngine::rcp_dynamic_cast<const SymEngine::DiagonalMatrix>(H0.second);
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
        // A density matrix is Hermitian and has a trace one
        auto val_rho0 = get_values(rho0.second);
        auto val_rho0_H = get_values(
            SymEngine::conjugate_matrix(SymEngine::transpose(rho0.second))
        );
        for (std::size_t i=0; i<val_rho0.size(); ++i) {
            //FIXME: SymEngine::neq may fail for comparing `conjugate`
            if (SymEngine::neq(*val_rho0[i], *SymEngine::simplify(val_rho0_H[i])))
                throw SymEngine::SymEngineException(
                    "Density matrix must be Hermitian: " + stringify(rho0.second)
                );
        }
        if (!is_zero_quantity(SymEngine::add(
            SymEngine::minus_one, SymEngine::trace(rho0.second)
        ))) throw SymEngine::SymEngineException(
                "Density matrix must have a trace one: " + stringify(rho0.second)
            );
        rho0_ = rho0;
        // Check the validity of field operators
        SymEngine::set_basic perturbations;
        for (const auto& oper: V_) {
            auto dependencies = oper.first->get_dependencies();
            if (dependencies.size()!=1) throw SymEngine::SymEngineException(
                "Each field operator should depend only on one perturbation: "
                + stringify(oper.first)
            );
            if (perturbations.find(dependencies.begin()->first)!=perturbations.end())
                throw SymEngine::SymEngineException(
                    "Field operators should depend on different perturbations: "
                    + stringify(oper.first)
                );
            perturbations.insert(dependencies.begin()->first);
            //Elements of V's can be symbols that we cannot check the hermicity
            //// All V's must be Hermitian
            //if (SymEngine::neq(
            //    *oper.second,
            //    *SymEngine::conjugate_matrix(SymEngine::transpose(oper.second))
            //)) throw SymEngine::SymEngineException(
            //        "Each field operator must be Hermitian: " + stringify(oper.second)
            //    );
        }
        // Store unperturbed density matrix
        rho_cached_[0] = DensityDerivative({
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
                auto max_order_cached = rho_cached_.size()-1;
                // Return when derivatives of the density matrix already computed
                if (derivatives.size()<=max_order_cached)
                    for (const auto& iter: rho_cached_[derivatives.size()])
                        if (SymEngine::unified_eq(derivatives, iter.first))
                            return iter.second;
                // Find matched lower order derivatives of the density matrix
                unsigned int lower_order = derivatives.size()<=max_order_cached
                                         ? derivatives.size()-1 : max_order_cached;
                auto lower_derivatives = SymEngine::multiset_basic(
                    derivatives.begin(),
                    std::next(derivatives.begin(), lower_order)
                );
                auto lower_found = false;
                SymEngine::RCP<const SymEngine::MatrixExpr> rho_matrix;
                for (; lower_order>0; --lower_order) {
                    auto last_perturbation = std::next(
                        lower_derivatives.begin(), lower_order-1
                    );
                    for (const auto& iter: rho_cached_[lower_order]) {
                        if (SymEngine::unified_eq(iter.first, lower_derivatives)) {
                            rho_matrix = iter.second;
                            lower_found = true;
                            break;
                        }
                    }
                    if (lower_found) break;
                    lower_derivatives.erase(last_perturbation);
                }
                if (!lower_found) rho_matrix = rho_cached_[0].begin()->second;
                // Loop over lower order derivatives
                for (; lower_order<derivatives.size(); ++lower_order) {
                    // Find the field operator according to its perturbation
                    auto field_pert = *(std::next(derivatives.begin(), lower_order));
                    SymEngine::RCP<const SymEngine::MatrixExpr> field_matrix;
                    auto field_found = true;
                    for (const auto& oper: V_) {
                        if (SymEngine::eq(
                            *field_pert,
                            *oper.first->get_dependencies().begin()->first
                        )) {
                            field_matrix = oper.second;
                            field_found = false;
                            break;
                        }
                    }
                    if (field_found) throw SymEngine::SymEngineException(
                        "Invalid perturbation for the external field "
                        + stringify(field_pert)
                    );
                    // Compute higher order derivatives of the density matrix
                    auto higher_order = lower_order+1;
                    auto higher_derivatives = SymEngine::multiset_basic(
                        derivatives.begin(),
                        std::next(derivatives.begin(), higher_order)
                    );
                    auto freq_sum = get_frequency_sum(higher_derivatives);
                    auto val_V_rho = get_values(SymEngine::matrix_mul({
                        field_matrix, rho_matrix
                    }));
                    auto val_rho_V = get_values(SymEngine::matrix_mul({
                        rho_matrix, field_matrix
                    }));
                    SymEngine::vec_basic val_rho;
                    for (std::size_t i=0; i<omega_.size(); ++i) {
                        val_rho.push_back(SymEngine::simplify(
                            SymEngine::div(
                                SymEngine::mul(
                                    SymEngine::integer(higher_order),
                                    SymEngine::sub(val_V_rho[i], val_rho_V[i])
                                ),
                                SymEngine::sub(freq_sum, omega_[i])
                            )
                        ));
                    }
                    // Update cached derivatives of density matrix by inserting
                    // the computed higher order derivatives of the density matrix
                    rho_matrix = SymEngine::immutable_dense_matrix(
                        dim_operator_, dim_operator_, val_rho
                    );
                    auto iter = rho_cached_.find(higher_order);
                    if (iter==rho_cached_.end()) {
                        rho_cached_.insert(std::pair<unsigned int, DensityDerivative>{
                            higher_order,
                            DensityDerivative({
                                std::make_pair(higher_derivatives, rho_matrix)
                            })
                        });
                    }
                    else {
                        iter->second.push_back(
                            std::make_pair(higher_derivatives, rho_matrix)
                        );
                    }
                }
                return rho_matrix;
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
