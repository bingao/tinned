/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of function evaluator.

   2023-12-03, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>

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
//#include "Tinned/CompositeFunction.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    template<typename FunctionType>
    class FunctionEvaluator: public SymEngine::BaseVisitor<FunctionEvaluator<FunctionType>>
    {
        protected:
            FunctionType result_;

            virtual FunctionType eval_nonel_function(const NonElecFunction& x) = 0;
            virtual FunctionType eval_2el_energy(const TwoElecEnergy& x) = 0;
            virtual FunctionType eval_xc_energy(const ExchCorrEnergy& x) = 0;
            virtual FunctionType eval_trace(
                const SymEngine::RCP<const SymEngine::Basic>& scalar,
                const SymEngine::vec_basic& factors
            ) = 0;
            virtual void eval_fun_addition(FunctionType& f, const FunctionType& g) = 0;
            virtual void eval_fun_scale(
                const SymEngine::RCP<const SymEngine::Number>& scalar,
                FunctionType& fun
            ) = 0;

        public:
            explicit FunctionEvaluator() = default;

            inline FunctionType apply(const SymEngine::RCP<const SymEngine::Basic>& x)
            {
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
                for (std::size_t i=1; i<args.size(); ++i)
                    eval_fun_addition(result_, apply(args[i]));
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
                    result_ = eval_nonel_function(
                        SymEngine::down_cast<const NonElecFunction&>(x)
                    );
                }
                else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
                    result_ = eval_2el_energy(
                        SymEngine::down_cast<const TwoElecEnergy&>(x)
                    );
                }
                //else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
                //    result_ = eval_composite_function(
                //        SymEngine::down_cast<const CompositeFunction&>(x)
                //    );
                //}
                else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
                    result_ = eval_xc_energy(
                        SymEngine::down_cast<const ExchCorrEnergy&>(x)
                    );
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
                auto arg = x.get_args()[0];
                if (SymEngine::is_a_sub<const SymEngine::MatrixMul>(*arg)) {
                    auto p = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(arg);
                    result_ = eval_trace(p->get_scalar(), p->get_factors());
                }
                else {
                    throw SymEngine::NotImplementedError(
                        "FunctionEvaluator::bvisit() not implemented for Trace "
                        + stringify(arg)
                    );
                }
            }
    };
}
