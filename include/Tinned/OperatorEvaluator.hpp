/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of operator evaluator.

   2024-06-14, Bin Gao:
   * add member variable `derivatives_` to hold derivatives of symbols to
     evaluate so that users can, for example, consider only interesting
     components of perturbations.

   2023-12-01, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>
#include <vector>

#include <symengine/basic.h>
#include <symengine/constants.h>
#include <symengine/number.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/transpose.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

#include "Tinned/PerturbedParameter.hpp"
#include "Tinned/ConjugateTranspose.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

//#include "Tinned/AdjointMap.hpp"
//#include "Tinned/ClusterConjHamiltonian.hpp"

#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    template<typename OperatorType>
    class OperatorEvaluator: public SymEngine::BaseVisitor<OperatorEvaluator<OperatorType>>
    {
        protected:
            // Each symbol adds its derivative to the end of the vector, and
            // its parent checks the validity of derivative and may remove it
            // when the child symbol has been evaluated.
            std::vector<SymEngine::multiset_basic> derivatives_;
            OperatorType result_;

            virtual OperatorType eval_pert_parameter(const PerturbedParameter& x) = 0;
            // return A^{\dagger}
            virtual OperatorType eval_hermitian_transpose(const OperatorType& A) = 0;
            virtual OperatorType eval_1el_density(const OneElecDensity& x) = 0;
            virtual OperatorType eval_1el_operator(const OneElecOperator& x) = 0;
            virtual OperatorType eval_2el_operator(const TwoElecOperator& x) = 0;
            virtual OperatorType eval_temporum_operator(const TemporumOperator& x) = 0;
            virtual OperatorType eval_temporum_overlap(const TemporumOverlap& x) = 0;
            // return A^{*}
            virtual OperatorType eval_conjugate_matrix(const OperatorType& A) = 0;
            // return A^{T}
            virtual OperatorType eval_transpose(const OperatorType& A) = 0;
            // A = A + B
            virtual void eval_oper_addition(OperatorType& A, const OperatorType& B) = 0;
            // return A * B
            virtual OperatorType eval_oper_multiplication(
                const OperatorType& A, const OperatorType& B
            ) = 0;
            // A = scalar * A
            virtual void eval_oper_scale(
                const SymEngine::RCP<const SymEngine::Number>& scalar,
                OperatorType& A
            ) = 0;

            // Update the derivative of a multiplication
            inline void update_mul_derivative() noexcept
            {
                auto derivative = derivatives_.back();
                derivatives_.pop_back();
                // The derivative of a multiplication is the union of
                // derivatives of all its factors
                derivatives_.back().insert(derivative.begin(), derivative.end());
            }

        public:
            explicit OperatorEvaluator() = default;

            inline OperatorType apply(const SymEngine::RCP<const SymEngine::Basic>& x)
            {
                derivatives_.clear();
                x->accept(*this);
                return result_;
            }

            inline std::vector<SymEngine::multiset_basic> get_derivatives() const
            {
                return derivatives_;
            }

            void bvisit(const SymEngine::Basic& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::bvisit() not implemented for " + stringify(x)
                );
            }

            void bvisit(const SymEngine::MatrixSymbol& x)
            {
                if (SymEngine::is_a_sub<const PerturbedParameter>(x)) {
                    auto& op = SymEngine::down_cast<const PerturbedParameter&>(x);
                    derivatives_.push_back(op.get_derivatives());
                    result_ = eval_pert_parameter(op);
                }
                else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
                    auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
                    result_ = eval_hermitian_transpose(apply(op.get_arg()));
                }
                else if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
                    auto& op = SymEngine::down_cast<const OneElecDensity&>(x);
                    derivatives_.push_back(op.get_derivatives());
                    result_ = eval_1el_density(op);
                }
                else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
                    auto& op = SymEngine::down_cast<const OneElecOperator&>(x);
                    derivatives_.push_back(op.get_derivatives());
                    result_ = eval_1el_operator(op);
                }
                else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
                    auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
                    auto op_derivatives = op.get_derivatives();
                    auto state_derivatives = op.get_state()->get_derivatives();
                    op_derivatives.insert(
                        state_derivatives.begin(), state_derivatives.end()
                    );
                    derivatives_.push_back(op_derivatives);
                    result_ = eval_2el_operator(op);
                }
                else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
                    auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
                    derivatives_.push_back(op.get_derivatives());
                    result_ = eval_xc_potential(op);
                }
                else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
                    auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
                    derivatives_.push_back(op.get_derivatives());
                    result_ = eval_temporum_operator(op);
                }
                else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
                    auto& op = SymEngine::down_cast<const TemporumOverlap&>(x);
                    derivatives_.push_back(op.get_derivatives());
                    result_ = eval_temporum_overlap(op);
                }
                else {
                    throw SymEngine::NotImplementedError(
                        "OperatorEvaluator::bvisit() not implemented for MatrixSymbol "
                        + stringify(x)
                    );
                }
            }

            void bvisit(const SymEngine::ConjugateMatrix& x)
            {
                result_ = eval_conjugate_matrix(apply(x.get_arg()));
            }

            void bvisit(const SymEngine::Transpose& x)
            {
                result_ = eval_transpose(apply(x.get_arg()));
            }

            void bvisit(const SymEngine::MatrixAdd& x)
            {
                auto args = x.get_args();
                result_ = apply(args[0]);
                for (std::size_t i=1; i<args.size(); ++i) {
                    auto val = apply(args[i]);
                    // Arguments of `MatrixAdd` should have the same derivative
                    if (SymEngine::unified_eq(
                        derivatives_.back(), derivatives_[derivatives_.size()-2]
                    )) {
                        eval_oper_addition(result_, val);
                        // We keep only the derivative of the first argument,
                        // which represents the derivative of `MatrixAdd`
                        derivatives_.pop_back();
                    }
                    else {
                        throw SymEngine::NotImplementedError(
                            "OperatorEvaluator::bvisit() got invalid MatrixAdd "
                            + stringify(x)
                        );
                    }
                }
            }

            void bvisit(const SymEngine::MatrixMul& x)
            {
                auto factors = x.get_factors();
                switch (factors.size()) {
                    case 1:
                        result_ = apply(factors[0]);
                        break;
                    case 2:
                        auto val0 = apply(factors[0]);
                        auto val1 = apply(factors[1]);
                        result_ = eval_oper_multiplication(val0, val1);
                        update_mul_derivative();
                        break;
                    case 3:
                        auto val0 = apply(factors[0]);
                        result_ = apply(factors[1]);
                        auto val1 = eval_oper_multiplication(val0, result_);
                        update_mul_derivative();
                        val0 = apply(factors[2]);
                        result_ = eval_oper_multiplication(val1, val0);
                        update_mul_derivative();
                        break;
                    default:
                        auto val0 = apply(factors[0]);
                        auto val = apply(factors[1]);
                        auto val1 = eval_oper_multiplication(val0, val);
                        update_mul_derivative();
                        for (std::size_t i=2; i<factors.size()-1; ++i) {
                            val = apply(factors[i]);
                            if (i%2==0) {
                                val0 = eval_oper_multiplication(val1, val);
                            }
                            else {
                                val1 = eval_oper_multiplication(val0, val);
                            }
                            update_mul_derivative();
                        }
                        val = apply(factors.back());
                        if (factors.size()%2==0) {
                            result_ = eval_oper_multiplication(val0, val);
                        }
                        else {
                            result_ = eval_oper_multiplication(val1, val);
                        }
                        update_mul_derivative();
                        break;
                }
                auto scalar = x.get_scalar();
                if (SymEngine::neq(*scalar, *SymEngine::one)) {
                    if (SymEngine::is_a_Number(*scalar)) {
                        eval_oper_scale(
                            SymEngine::rcp_dynamic_cast<const SymEngine::Number>(scalar),
                            result_
                        );
                    }
                    else {
                        throw SymEngine::NotImplementedError(
                            "OperatorEvaluator::bvisit() not implemented for scalar "
                            + stringify(scalar)
                        );
                    }
                }
            }
    };
}
