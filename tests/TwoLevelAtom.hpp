/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of two-level atom system.

   2024-07-15, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/number.h>
#include <symengine/integer.h>
#include <symengine/constants.h>
#include <symengine/derivative.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/matrices/immutable_dense_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"
#include "Tinned/OperatorEvaluator.hpp"
#include "Tinned/FunctionEvaluator.hpp"

// This two-level atom system is adopted from section 5.10.3, "Principles and
// Practices of Molecular Properties: Theory, Modeling and Simulations",
// Patrick Norman, Kenneth Ruud and Trond Saue.
namespace Tinned
{
    // Evaluator for different (electron) operators
    class TwoLevelOperator: public OperatorEvaluator<SymEngine::RCP<const SymEngine::MatrixExpr>>
    {
        protected:
            // Maximum allowed order
            unsigned int max_order_;
            // Unperturbed Hamiltonian and its value
            std::pair<SymEngine::RCP<const OneElecOperator>,
                      SymEngine::RCP<const SymEngine::MatrixExpr>> H0_;
            // Density matrix and the value of unperturbed one
            std::pair<SymEngine::RCP<const OneElecOperator>,
                      SymEngine::RCP<const SymEngine::MatrixExpr>> rho_;
            // External field's operators and their values, each operator
            // should depend only on one unique perturbation
            std::map<SymEngine::RCP<const OneElecOperator>,
                     SymEngine::RCP<const SymEngine::MatrixExpr>,
                     SymEngine::RCPBasicKeyLess> V_;
            // Transition angular frequency matrix
            SymEngine::RCP<const SymEngine::ImmutableDenseMatrix> omega_;
            // Cached derivatives of density matrix <order, <perturbations, derivatives>>
            //FIXME: not optimal for identical perturbations
            typedef std::vector<std::pair<SymEngine::vec_basic, SymEngine::RCP<const SymEngine::MatrixExpr>>>
                DensityDerivative;
            std::map<unsigned int, DensityDerivative> rho_all_derivatives_;

            SymEngine::vec_basic eval_pert_parameter(const PerturbedParameter& x) override
            {
                throw SymEngine::SymEngineException(
                    "PerturbedParameter is not allowed: " + stringify(x)
                );
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_hermitian_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override
            {
                auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
                SymEngine::vec_basic values;
                for (std::size_t j=0; j<op->ncols(); ++j)
                    for (std::size_t i=0; i<op->nrows(); ++i)
                        values.push_back(op->get(i, j)->conjugate());
                return SymEngine::immutable_dense_matrix(op->ncols(), op->nrows(), values);
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_1el_density(const OneElecDensity& x) override
            {
                auto derivatives = x.get_derivatives();
                if (derivatives.size()>max_order_) SymEngine::SymEngineException(
                    "Invalid order: " + std::to_string(derivatives.size())
                    + "(" + std::to_string(max_order_) + ")"
                );
                switch (derivatives.size()) {
                    case 0:
                        return rho_.second;
                    default:
                        if (derivatives.size()>rho_all_derivatives_.size()) {
                            for (unsigned int order=rho_all_derivatives_.size();
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
                                    // external field, and we need to find the
                                    // value of its corresponding operator
                                    SymEngine::RCP<const SymEngine::MatrixExpr> oper_field;
                                    bool not_found = true;
                                    for (const auto& oper: V_) {
                                        if (SymEngine::eq(
                                            permut_derivatives.back(),
                                            oper.first->get_dependencies().begin()->first
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
                                    not_found = true;
                                    // Try to find the matched lower order
                                    // derivatives of density matrix
                                    for (const auto& lower_derivatives: rho_all_derivatives_[order]) {
                                        bool pert_matched = true;
                                        for (std::size_t i=0; i<lower_derivatives.first.size(); ++i)
                                            if (SymEngine::neq(
                                                permut_derivatives[i], lower_derivatives.first[i]
                                            )) {
                                                pert_matched = false;
                                                break;
                                            }
                                        if (pert_matched) {
                                            auto rho_commutator
                                                = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
                                                      SymEngine::matrix_add({
                                                          SymEngine::matrix_mul({
                                                              oper_field,
                                                              lower_derivatives.second
                                                          }),
                                                          SymEngine::matrix_mul({
                                                              SymEngine::minus_one,
                                                              lower_derivatives.second,
                                                              oper_field
                                                          })
                                                      })
                                                  );
                                            SymEngine::vec_basic val_rho;
                                            for (std::size_t i=0; i<omega_->nrows(); ++i)
                                                for (std::size_t j=0; j<omega_->ncols(); ++j)
                                                    val_rho.push_back(
                                                        SymEngine::div(
                                                            rho_commutator->get(i, j),
                                                            SymEngine::sub(
                                                                freq_sum,
                                                                omega_->get(i, j)
                                                            )
                                                        )
                                                    );
                                            rho_derivatives.push_back(
                                                std::make_pair(
                                                    permut_derivatives,
                                                    SymEngine::immutable_dense_matrix(
                                                        omega_->nrows(),
                                                        omega_->ncols(),
                                                        val_rho
                                                    )
                                                )
                                            );
                                            not_found = false;
                                            break;
                                        }
                                    }
                                    if (not_found) throw SymEngine::SymEngineException(
                                        "Invalid ower order derivatives of density matrix "
                                        + stringify(permut_derivatives)
                                    );
                                } while (remaining);
                                rho_all_derivatives_[order+1] = rho_derivatives;
                            }
                        }
                        for (const auto& rho_derivatives: rho_all_derivatives_[derivatives.size()])
                            if (SymEngine::unified_eq(derivatives, rho_derivatives.first))
                                return rho_derivatives.second;
                        throw SymEngine::SymEngineException(
                            "Density matrix with invalid derivatives : " + stringify(x)
                        );
                }
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_1el_operator(const OneElecOperator& x) override
            {
                if (x.get_name())
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_2el_operator(const TwoElecOperator& x) override
            {
                throw SymEngine::SymEngineException(
                    "TwoElecOperator is not allowed: " + stringify(x)
                );
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_temporum_operator(const TemporumOperator& x) override
            {
                result_ = apply(x.get_target());
                eval_oper_scale(x.get_frequency(), result_);
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_temporum_overlap(const TemporumOverlap& x) override
            {
                throw SymEngine::SymEngineException(
                    "TemporumOverlap is not allowed: " + stringify(x)
                );
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_conjugate_matrix(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override
            {
                auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
                SymEngine::vec_basic values;
                for (std::size_t i=0; i<op->nrows(); ++i)
                    for (std::size_t j=0; j<op->ncols(); ++j)
                        values.push_back(op->get(i, j)->conjugate());
                return SymEngine::immutable_dense_matrix(op->nrows(), op->ncols(), values);
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override
            {
                auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
                SymEngine::vec_basic values;
                for (std::size_t j=0; j<op->ncols(); ++j)
                    for (std::size_t i=0; i<op->nrows(); ++i)
                        values.push_back(op->get(i, j));
                return SymEngine::immutable_dense_matrix(op->ncols(), op->nrows(), values);
            }

            void eval_oper_addition(
                SymEngine::RCP<const SymEngine::MatrixExpr>& A,
                const SymEngine::RCP<const SymEngine::MatrixExpr>& B
            ) override
            {
                auto op_A = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
                auto op_B = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(B);
                SYMENGINE_ASSERT(
                    op_A->nrows()==op_B->nrows() && op_A->ncols()==op_B->ncols()
                )
                SymEngine::vec_basic values;
                for (std::size_t i=0; i<op_A->nrows(); ++i)
                    for (std::size_t j=0; j<op_A->ncols(); ++j)
                        values.push_back(
                            SymEngine::add(op_A->get(i, j), op_B->get(i, j))
                        );
                A = SymEngine::immutable_dense_matrix(op_A->nrows(), op_A->ncols(), values);
            }

            SymEngine::RCP<const SymEngine::MatrixExpr> eval_oper_multiplication(
                const SymEngine::RCP<const SymEngine::MatrixExpr>& A,
                const SymEngine::RCP<const SymEngine::MatrixExpr>& B
            ) override
            {
                return SymEngine::matrix_mul({A, B});
            }

            void eval_oper_scale(
                const SymEngine::RCP<const SymEngine::Number>& scalar,
                SymEngine::RCP<const SymEngine::MatrixExpr>& A
            ) override
            {
                auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
                SymEngine::vec_basic values;
                for (std::size_t i=0; i<op->nrows(); ++i)
                    for (std::size_t j=0; j<op->ncols(); ++j)
                        values.push_back(SymEngine::mul(scalar, op->get(i, j)));
                A = SymEngine::immutable_dense_matrix(op->nrows(), op->ncols(), values);
            }

        public:
            explicit TwoLevelOperator(
                const SymEngine::RCP<const SymEngine::MatrixExpr>& H0,
                const std::map<SymEngine::RCP<const Perturbation>,
                               SymEngine::RCP<const SymEngine::MatrixExpr>,
                               SymEngine::RCPBasicKeyLess()>& V_mapping,
                const SymEngine::RCP<const SymEngine::MatrixExpr>& rho0
            ) : max_order_(7), H0_(H0), V_mapping_(V_mapping)
            {
                // Density matrix should be idempotent and has purity one
                auto Z = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
                    SymEngine::matrix_add({
                        SymEngine::matrix_mul({rho0, rho0}),
                        SymEngine::matrix_mul({SymEngine::minus_one, rho0})
                    })
                );
                for (const auto& val: Z->get_values())
                    if (!is_zero_quantity(val)) throw SymEngine::SymEngineException(
                        "Density matrix isn't idempotent: " + stringify(rho0)
                    );
                if (!is_zero_quantity(SymEngine::add(SymEngine::minus_one, SymEngine::trace(rho0))))
                    throw SymEngine::SymEngineException(
                        "Density matrix doesn't have purity one: " + stringify(rho0)
                    );
                rho0_ = rho0;
            }

            ~TwoLevelOperator() = default;
    }

    // Evaluator for different expectation values
    class TwoLevelFunction: public FunctionEvaluator<SymEngine::RCP<const SymEngine::Basic>>
    {
        private:
            SymEngine::vec_basic perturbations_;

        protected:
            SymEngine::RCP<const SymEngine::Basic> eval_nonel_function(
                const NonElecFunction& x
            ) override
            {
            }

            SymEngine::RCP<const SymEngine::Basic> eval_2el_energy(
                const TwoElecEnergy& x
            ) override
            {
            }

            SymEngine::RCP<const SymEngine::Basic> eval_xc_energy(
                const ExchCorrEnergy& x
            ) override
            {
            }

            SymEngine::RCP<const SymEngine::Basic> eval_trace(
                const SymEngine::RCP<const SymEngine::Basic>& scalar,
                const SymEngine::vec_basic& factors
            ) override
            {
            }

            void eval_fun_addition(
                SymEngine::RCP<const SymEngine::Basic>& f,
                const SymEngine::RCP<const SymEngine::Basic>& g
            ) override
            {
                f = SymEngine::add(f, g);
            }

            void eval_fun_scale(
                const SymEngine::RCP<const SymEngine::Number>& scalar,
                SymEngine::RCP<const SymEngine::Basic>& fun
            ) override
            {
                fun = SymEngine::mul(scalar, fun);
            }

        public:
            explicit TwoLevelFunction(const SymEngine::vec_basic& perturbations):
                perturbations_(perturbations) {}
            ~TwoLevelFunction() = default;
    };
}
