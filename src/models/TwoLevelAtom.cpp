#include <cstddef>
#include <iterator>
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
#include <symengine/simplify.h>

#include "Tinned/models/TwoLevelAtom.hpp"
#include "Tinned/visitors/StringifyVisitor.hpp"

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
    TwoLevelOperator::eval_1el_density(const SymEngine::multiset_basic& derivatives)
    {
        // Higher order derivatives of the density matrix
        SymEngine::vec_basic val_rho;
        for (std::size_t i=0; i<dim_operator_*dim_operator_; ++i)
            val_rho.push_back(SymEngine::zero);
        // Loop over perturbation of the field operator
        auto field_pert=derivatives.begin();
        for (; field_pert!=derivatives.end(); ++field_pert) {
            // Find the field operator according to its perturbation
            SymEngine::RCP<const SymEngine::MatrixExpr> field_matrix;
            auto field_found = true;
            for (const auto& oper: V_) {
                if (SymEngine::eq(
                    *(*field_pert),
                    *oper.first->get_dependencies().begin()->first
                )) {
                    field_matrix = oper.second;
                    field_found = false;
                    break;
                }
            }
            if (field_found) throw SymEngine::SymEngineException(
                "Invalid perturbation for the external field " + stringify(*field_pert)
            );
            // Lower order derivatives of the density matrix
            auto lower_derivatives = SymEngine::multiset_basic(
                derivatives.begin(), field_pert
            );
            lower_derivatives.insert(std::next(field_pert, 1), derivatives.end());
            SymEngine::RCP<const SymEngine::MatrixExpr> rho_lower_matrix;
            auto lower_order = derivatives.size()-1;
            auto rho_lower = rho_cached_.find(lower_order);
            if (rho_lower==rho_cached_.end()) {
                rho_lower_matrix = eval_1el_density(lower_derivatives);
            }
            else {
                auto not_found = true;
                for (const auto& iter: rho_lower->second) {
                    if (SymEngine::unified_eq(iter.first, lower_derivatives)) {
                        rho_lower_matrix = iter.second;
                        not_found = false;
                        break;
                    }
                }
                if (not_found) rho_lower_matrix = eval_1el_density(lower_derivatives);
            }
            // Commutator of the field operator and lower order derivatives of
            // density matrix
            auto val_V_rho = get_values(SymEngine::matrix_mul({
                field_matrix, rho_lower_matrix
            }));
            auto val_rho_V = get_values(SymEngine::matrix_mul({
                rho_lower_matrix, field_matrix
            }));
            for (std::size_t i=0; i<val_rho.size(); ++i)
                val_rho[i] = SymEngine::simplify(SymEngine::add(
                    val_rho[i], SymEngine::sub(val_V_rho[i], val_rho_V[i])
                ));
        }
        // Compute higher order derivatives of the density matrix
        auto freq_sum = get_frequency_sum(derivatives);
        for (std::size_t i=0; i<val_rho.size(); ++i)
            val_rho[i] = SymEngine::div(val_rho[i], SymEngine::sub(freq_sum, omega_[i]));
        auto rho_matrix = SymEngine::immutable_dense_matrix(
            dim_operator_, dim_operator_, val_rho
        );
        // Update cached derivatives of density matrix by inserting the
        // computed higher order derivatives of the density matrix
        auto rho_higher = rho_cached_.find(derivatives.size());
        if (rho_higher==rho_cached_.end()) {
            rho_cached_.insert(std::pair<unsigned int, DensityDerivative>{
                derivatives.size(),
                DensityDerivative({std::make_pair(derivatives, rho_matrix)})
            });
        }
        else {
            rho_higher->second.push_back(std::make_pair(derivatives, rho_matrix));
        }
        return rho_matrix;
    }

    //SymEngine::RCP<const SymEngine::MatrixExpr>
    //TwoLevelOperator::eval_hermitian_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A)
    //{
    //    auto conj = SymEngine::conjugate_matrix(A);
    //    return SymEngine::transpose(conj);
    //}

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
                return eval_1el_density(derivatives);
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
                    if (derivatives.size()==1 &&
                        SymEngine::eq(
                            *(*derivatives.begin()),
                            *(oper.first->get_dependencies().begin()->first)
                        )
                    ) {
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

    //SymEngine::RCP<const SymEngine::MatrixExpr>
    //TwoLevelOperator::eval_temporum_operator(const TemporumOperator& x)
    //{
    //    result_ = apply(x.get_target());
    //    eval_oper_scale(x.get_frequency(), result_);
    //    return result_;
    //}
    //
    //SymEngine::RCP<const SymEngine::MatrixExpr> TwoLevelOperator::eval_conjugate_matrix(
    //    const SymEngine::RCP<const SymEngine::MatrixExpr>& A
    //)
    //{
    //    return SymEngine::conjugate_matrix(A);
    //}

    //SymEngine::RCP<const SymEngine::MatrixExpr>
    //TwoLevelOperator::eval_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A)
    //{
    //    return SymEngine::transpose(A);
    //}

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
