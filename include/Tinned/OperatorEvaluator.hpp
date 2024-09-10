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

            virtual OperatorType eval_pert_parameter(const PerturbedParameter& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_pert_parameter() is not implemented"
                );
            }

            // return A^{\dagger}
            virtual OperatorType eval_hermitian_transpose(const OperatorType& A)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_hermitian_transpose() is not implemented"
                );
            }

            virtual OperatorType eval_1el_density(const OneElecDensity& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_1el_density() is not implemented"
                );
            }

            virtual OperatorType eval_1el_operator(const OneElecOperator& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_1el_operator() is not implemented"
                );
            }

            virtual OperatorType eval_2el_operator(const TwoElecOperator& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_2el_operator() is not implemented"
                );
            }

            virtual OperatorType eval_xc_potential(const ExchCorrPotential& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_xc_potential() is not implemented"
                );
            }

            virtual OperatorType eval_temporum_operator(const TemporumOperator& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_temporum_operator() is not implemented"
                );
            }

            virtual OperatorType eval_temporum_overlap(const TemporumOverlap& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_temporum_overlap() is not implemented"
                );
            }

            // return A^{*}
            virtual OperatorType eval_conjugate_matrix(const OperatorType& A)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_conjugate_matrix() is not implemented"
                );
            }

            // return A^{T}
            virtual OperatorType eval_transpose(const OperatorType& A)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_transpose() is not implemented"
                );
            }

            // A = A + B
            virtual void eval_oper_addition(OperatorType& A, const OperatorType& B)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_oper_addition() is not implemented"
                );
            }

            // return A * B
            virtual OperatorType eval_oper_multiplication(
                const OperatorType& A, const OperatorType& B
            )
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_oper_multiplication() is not implemented"
                );
            }

            // A = scalar * A
            virtual void eval_oper_scale(
                const SymEngine::RCP<const SymEngine::Basic>& scalar,
                OperatorType& A
            )
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::eval_oper_scale() is not implemented"
                );
            }

            // Update the derivative of a multiplication of two factors
            inline void update_mul_derivative() noexcept
            {
                auto derivative = derivatives_.back();
                derivatives_.pop_back();
                // The derivative of a multiplication is the union of
                // derivatives of its factors
                if (derivative.size()>0) derivatives_.back().insert(
                    derivative.begin(), derivative.end()
                );
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
                    op.get_arg()->accept(*this);
                    result_ = eval_hermitian_transpose(result_);
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
                x.get_arg()->accept(*this);
                result_ = eval_conjugate_matrix(result_);
            }

            void bvisit(const SymEngine::Transpose& x)
            {
                x.get_arg()->accept(*this);
                result_ = eval_transpose(result_);
            }

            void bvisit(const SymEngine::MatrixAdd& x)
            {
                auto args = x.get_args();
                args[0]->accept(*this);
                auto sum = result_;
                for (std::size_t i=1; i<args.size(); ++i) {
                    args[i]->accept(*this);
                    // Arguments of `MatrixAdd` should have the same derivative
                    if (SymEngine::unified_eq(
                        derivatives_.back(), derivatives_[derivatives_.size()-2]
                    )) {
                        eval_oper_addition(sum, result_);
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
                result_ = sum;
            }

            void bvisit(const SymEngine::MatrixMul& x)
            {
                auto factors = x.get_factors();
                switch (factors.size()) {
                    case 1: {
                        factors[0]->accept(*this);
                        break;
                    }
                    case 2: {
                        factors[0]->accept(*this);
                        auto val0 = result_;
                        factors[1]->accept(*this);
                        auto val1 = result_;
                        result_ = eval_oper_multiplication(val0, val1);
                        update_mul_derivative();
                        break;
                    }
                    case 3: {
                        factors[0]->accept(*this);
                        auto val0 = result_;
                        factors[1]->accept(*this);
                        auto val1 = eval_oper_multiplication(val0, result_);
                        update_mul_derivative();
                        factors[2]->accept(*this);
                        val0 = result_;
                        result_ = eval_oper_multiplication(val1, val0);
                        update_mul_derivative();
                        break;
                    }
                    default: {
                        factors[0]->accept(*this);
                        auto val0 = result_;
                        factors[1]->accept(*this);
                        auto val1 = eval_oper_multiplication(val0, result_);
                        update_mul_derivative();
                        for (std::size_t i=2; i<factors.size()-1; ++i) {
                            factors[i]->accept(*this);
                            if (i%2==0) {
                                val0 = eval_oper_multiplication(val1, result_);
                            }
                            else {
                                val1 = eval_oper_multiplication(val0, result_);
                            }
                            update_mul_derivative();
                        }
                        factors.back()->accept(*this);
                        if (factors.size()%2==0) {
                            val1 = result_;
                            result_ = eval_oper_multiplication(val0, val1);
                        }
                        else {
                            val0 = result_;
                            result_ = eval_oper_multiplication(val1, val0);
                        }
                        update_mul_derivative();
                        break;
                    }
                }
                auto scalar = x.get_scalar();
                if (SymEngine::neq(*scalar, *SymEngine::one))
                    eval_oper_scale(scalar, result_);
            }

            virtual ~OperatorEvaluator() = default;
    };
}
