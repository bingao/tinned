/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of operator evaluator.

   2023-12-01, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>

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

#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"
#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    template<typename OperatorType>
    class OperatorEvaluator: public SymEngine::BaseVisitor<OperatorEvaluator<OperatorType>>
    {
        protected:
            bool in_place_conjugate_;
            bool in_place_transpose_;
            OperatorType result_;

            virtual OperatorType eval_1el_density(const OneElecDensity& x) = 0;
            virtual OperatorType eval_1el_operator(const OneElecOperator& x) = 0;
            virtual OperatorType eval_2el_operator(const TwoElecOperator& x) = 0;
            virtual OperatorType eval_temporum_operator(const TemporumOperator& x) = 0;
            virtual OperatorType eval_temporum_overlap(const TemporumOverlap& x) = 0;
            virtual void eval_conjugate_matrix(OperatorType& A, OperatorType& B) = 0;
            virtual void eval_transpose(OperatorType& A, OperatorType& B) = 0;
            virtual void eval_oper_addition(OperatorType& A, const OperatorType& B) = 0;
            virtual void eval_oper_multiplication(
                OperatorType& A, const OperatorType& B
            ) = 0;
            virtual void eval_oper_scale(
                const SymEngine::RCP<const SymEngine::Number>& scalar,
                OperatorType& A
            ) = 0;

        public:
            explicit OperatorEvaluator(
                const bool inPlaceConjugate = false,
                const bool inPlaceTranspose = false
            ): in_place_conjugate_(inPlaceConjugate),
               in_place_transpose_(inPlaceTranspose) {}

            inline OperatorType apply(const SymEngine::RCP<const SymEngine::Basic>& x)
            {
                x->accept(*this);
                return result_;
            }

            void bvisit(const SymEngine::Basic& x)
            {
                throw SymEngine::NotImplementedError(
                    "OperatorEvaluator::bvisit() not implemented for " + stringify(x)
                );
            }

            void bvisit(const SymEngine::MatrixSymbol& x)
            {
                if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
                    result_ = eval_1el_density(
                        SymEngine::down_cast<const OneElecDensity&>(x)
                    );
                }
                else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
                    result_ = eval_1el_operator(
                        SymEngine::down_cast<const OneElecOperator&>(x)
                    );
                }
                else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
                    result_ = eval_2el_operator(
                        SymEngine::down_cast<const TwoElecOperator&>(x)
                    );
                }
                else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
                    result_ = eval_xc_potential(
                        SymEngine::down_cast<const ExchCorrPotential&>(x)
                    );
                }
                else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
                    result_ = eval_temporum_operator(
                        SymEngine::down_cast<const TemporumOperator&>(x)
                    );
                }
                else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
                    result_ = eval_temporum_overlap(
                        SymEngine::down_cast<const TemporumOverlap&>(x)
                    );
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
                if (in_place_conjugate_) {
                    result_ = apply(x.get_arg());
                    eval_conjugate_matrix(result_, result_);
                }
                else {
                    eval_conjugate_matrix(apply(x.get_arg()), result_);
                }
            }

            void bvisit(const SymEngine::Transpose& x)
            {
                if (in_place_transpose_) {
                    result_ = apply(x.get_arg());
                    eval_transpose(result_, result_);
                }
                else {
                    eval_transpose(apply(x.get_arg()), result_);
                }
            }

            void bvisit(const SymEngine::MatrixAdd& x)
            {
                auto args = x.get_args();
                result_ = apply(args[0]);
                for (std::size_t i=1; i<args.size(); ++i)
                    eval_oper_addition(result_, apply(args[i]));
            }

            void bvisit(const SymEngine::MatrixMul& x)
            {
                auto factors = x.get_factors();
                result_ = apply(factors[0]);
                for (std::size_t i=1; i<factors.size(); ++i)
                    eval_oper_multiplication(result_, apply(factors[i]));
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
