/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of cleaning `TemporumOperator` objects.

   2024-05-08, Bin Gao:
   * first version
*/

#pragma once

#include <functional>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/mul.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

#include "Tinned/ZeroOperator.hpp"
#include "Tinned/ZerosRemover.hpp"

namespace Tinned
{
    // `TemporumCleaner` replaces `TemporumOperator` objects with their targets
    // multiplied by sums of perturbation frequencies. Undifferentiated targets
    // and targets with zero sums, and undifferentiated `TemporumOverlap`
    // objects will be set as zero quantities.
    class TemporumCleaner: public SymEngine::BaseVisitor<TemporumCleaner>
    {
        protected:
            SymEngine::RCP<const SymEngine::Basic> result_;

            // Template method for one argument function like classes
            template<typename Fun, typename Arg>
            inline void clean_one_arg_f(
                Fun& x,
                const SymEngine::RCP<Arg>& arg,
                std::function<SymEngine::RCP<const SymEngine::Basic>(
                    const SymEngine::RCP<Arg>&
                )> constructor,
                const bool is_operator = true
            )
            {
                // We check only if the argument is/has `TemporumOperator` object
                auto new_arg = apply(arg);
                if (is_zero_quantity(*new_arg)) {
                    if (is_operator) {
                        result_ = make_zero_operator();
                    }
                    else {
                        result_ = SymEngine::zero;
                    }
                }
                else {
                    if (SymEngine::eq(*arg, *new_arg)) {
                        result_ = x.rcp_from_this();
                    }
                    else {
                        result_ = constructor(SymEngine::rcp_dynamic_cast<Arg>(new_arg));
                    }
                }
            }

        public:
            explicit TemporumCleaner() noexcept = default;

            inline SymEngine::RCP<const SymEngine::Basic> apply(
                const SymEngine::RCP<const SymEngine::Basic>& x
            )
            {
                x->accept(*this);
                return result_;
            }

            void bvisit(const SymEngine::Basic& x);
            void bvisit(const SymEngine::Add& x);
            void bvisit(const SymEngine::Mul& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            void bvisit(const SymEngine::Trace& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixAdd& x);
            void bvisit(const SymEngine::MatrixMul& x);
            void bvisit(const SymEngine::MatrixDerivative& x);
    };

    // Helper function to clean `TemporumOperator` and unperturbed
    // `TemporumOverlap` objects in `x`
    inline SymEngine::RCP<const SymEngine::Basic> clean_temporum(
        const SymEngine::RCP<const SymEngine::Basic>& x
    )
    {
        TemporumCleaner visitor;
        // Remove zero quantities
        auto result = remove_zeros(visitor.apply(x));
        if (result.is_null()) {
            if (SymEngine::is_a_sub<const SymEngine::MatrixExpr>(*x)) {
                return make_zero_operator();
            }
            else {
                return SymEngine::zero;
            }
        }
        else {
            return result;
        }
    }
}
