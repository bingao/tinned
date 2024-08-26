/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of function evaluator.

   2024-06-14, Bin Gao:
   * add member variable `derivatives_` to hold derivatives of symbols to
     evaluate so that users can, for example, consider only interesting
     components of perturbations.

   2023-12-03, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/mul.h>
#include <symengine/number.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

#include "Tinned/TwoElecEnergy.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/StringifyVisitor.hpp"
#include "Tinned/OperatorEvaluator.hpp"

namespace Tinned
{
    template<typename FunctionType, typename OperatorType>
    class FunctionEvaluator: public SymEngine::BaseVisitor<FunctionEvaluator<FunctionType, OperatorType>>
    {
        protected:
            // Each symbol adds its derivative to the end of the vector, and
            // its parent checks the validity of derivative and may remove it
            // when the child symbol has been evaluated.
            std::vector<SymEngine::multiset_basic> derivatives_;
            FunctionType result_;
            std::shared_ptr<OperatorEvaluator<OperatorType>> oper_evaluator_;

            virtual FunctionType eval_nonel_function(const NonElecFunction& x)
            {
                throw SymEngine::NotImplementedError(
                    "FunctionEvaluator::eval_nonel_function() is not implemented"
                );
            }

            virtual FunctionType eval_2el_energy(const TwoElecEnergy& x)
            {
                throw SymEngine::NotImplementedError(
                    "FunctionEvaluator::eval_2el_energy() is not implemented"
                );
            }

            virtual FunctionType eval_xc_energy(const ExchCorrEnergy& x)
            {
                throw SymEngine::NotImplementedError(
                    "FunctionEvaluator::eval_xc_energy() is not implemented"
                );
            }

            // return the trace of `A`
            virtual FunctionType eval_trace(const OperatorType& A)
            {
                throw SymEngine::NotImplementedError(
                    "FunctionEvaluator::eval_trace() is not implemented"
                );
            }

            // `f` = `f` + `g`
            virtual void eval_fun_addition(FunctionType& f, const FunctionType& g)
            {
                throw SymEngine::NotImplementedError(
                    "FunctionEvaluator::eval_fun_addition() is not implemented"
                );
            }

            // `f` = `scalar` * `f`
            virtual void eval_fun_scale(
                const SymEngine::RCP<const SymEngine::Number>& scalar,
                FunctionType& f
            )
            {
                throw SymEngine::NotImplementedError(
                    "FunctionEvaluator::eval_fun_scale() is not implemented"
                );
            }

        public:
            explicit FunctionEvaluator(
                const std::shared_ptr<OperatorEvaluator<OperatorType>>& operEvaluator
            ): oper_evaluator_(operEvaluator) {}

            inline FunctionType apply(const SymEngine::RCP<const SymEngine::Basic>& x)
            {
                derivatives_.clear();
                x->accept(*this);
                return result_;
            }

            void bvisit(const SymEngine::Basic& x)
            {
                throw SymEngine::NotImplementedError(
                    "FunctionEvaluator::bvisit() not implemented for " + stringify(x)
                );
            }

            void bvisit(const SymEngine::Add& x)
            {
                auto args = x.get_args();
                result_ = apply(args[0]);
                for (std::size_t i=1; i<args.size(); ++i) {
                    auto val = apply(args[i]);
                    // Arguments of `Add` should have the same derivative
                    if (SymEngine::unified_eq(
                        derivatives_.back(), derivatives_[derivatives_.size()-2]
                    )) {
                        eval_fun_addition(result_, val);
                        // We keep only the derivative of the first argument,
                        // which represents the derivative of `Add`
                        derivatives_.pop_back();
                    }
                    else {
                        throw SymEngine::NotImplementedError(
                            "FunctionEvaluator::bvisit() got invalid Add "
                            + stringify(x)
                        );
                    }
                }
            }

            void bvisit(const SymEngine::Mul& x)
            {
                SymEngine::RCP<const SymEngine::Number> scalar = SymEngine::one;
                unsigned int num_non_numbers = 0;
                for (auto const& arg: x.get_args()) {
                    if (SymEngine::is_a_Number(*arg)) {
                        scalar = SymEngine::mulnum(
                            scalar,
                            SymEngine::rcp_dynamic_cast<const SymEngine::Number>(arg)
                        );
                    }
                    else {
                        ++num_non_numbers;
                        if (num_non_numbers>1) {
                            throw SymEngine::NotImplementedError(
                                "FunctionEvaluator::bvisit() not implemented for the argument "
                                + stringify(arg)
                                + " of the multiplication "
                                + stringify(x)
                            );
                        }
                        else {
                            result_ = apply(arg);
                        }
                    }
                }
                if (SymEngine::neq(*scalar, *SymEngine::one))
                    eval_fun_scale(scalar, result_);
            }

            void bvisit(const SymEngine::FunctionSymbol& x)
            {
                if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
                    auto& op = SymEngine::down_cast<const NonElecFunction&>(x);
                    derivatives_.push_back(op.get_derivatives());
                    result_ = eval_nonel_function(op);
                }
                else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
                    auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
                    auto op_derivatives = op.get_derivatives();
                    auto inner_derivatives = op.get_inner_state()->get_derivatives();
                    op_derivatives.insert(
                        inner_derivatives.begin(), inner_derivatives.end()
                    );
                    auto outer_derivatives = op.get_outer_state()->get_derivatives();
                    op_derivatives.insert(
                        outer_derivatives.begin(), outer_derivatives.end()
                    );
                    derivatives_.push_back(op_derivatives);
                    result_ = eval_2el_energy(op);
                }
                else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
                    auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
                    derivatives_.push_back(op.get_derivatives());
                    result_ = eval_xc_energy(op);
                }
                else {
                    throw SymEngine::NotImplementedError(
                        "FunctionEvaluator::bvisit() not implemented for FunctionSymbol "
                        + stringify(x)
                    );
                }
            }

            void bvisit(const SymEngine::Trace& x)
            {
                auto arg = oper_evaluator_->apply(x.get_args()[0]);
                auto oper_derivatives = oper_evaluator_->get_derivatives();
                derivatives_.insert(
                    derivatives_.end(), oper_derivatives.begin(), oper_derivatives.end()
                );
                result_ = eval_trace(arg);
            }
    };
}
