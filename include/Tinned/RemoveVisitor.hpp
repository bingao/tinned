/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of removing specific symbols.

   2023-10-19, Bin Gao:
   * class RemoveVisitor will remove matching symbols, while class KeepVisitor
     (moved into another file) will keep matching symbols while removing others

   2023-10-16, Bin Gao:
   * formal rules for the removal procedure

   2023-09-25, Bin Gao:
   * first version
*/

#pragma once

#include <functional>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/complex.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/rational.h>
#include <symengine/symbol.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

namespace Tinned
{
    // Removing symbols if they match any given ones
    //
    // Rules for the removal procedure are that,
    //
    // (1) Nothing left after removal, which is different from replacing
    //     symbols by zero -- but we may find out the effect of removal by
    //     setting symbols to remove as zero.
    //
    // (2) If a symbol is involved in an operation, the operation should still
    //     be valid/meaningful after removal. Otherwise, either the whole
    //     operation will be removed or such a removal is not allowed.
    //
    //     For example, a symbol removed from an addition or a multiplication
    //     is allowed and is the same as replacing it by zero. That is, the
    //     multiplication will be removed if one of its factor is removed.
    //
    //     Removing an exponent from a power, removing a variable that a
    //     derivative is with respect to, and removing an element from a
    //     matrix, all are not allowed, because the left power and
    //     differentiation operations, and the matrix become invalid.
    //
    // (3) Symbols and their derivatives are different for the removal
    //     procedure.
    class RemoveVisitor: public SymEngine::BaseVisitor<RemoveVisitor>
    {
        protected:
            SymEngine::RCP<const SymEngine::Basic> result_;
            const SymEngine::vec_basic& symbols_;
            std::function<bool(const SymEngine::Basic&)> condition_;

            // Check equality for `x` and symbols to be removed
            inline bool is_equal(const SymEngine::Basic& x) {
                for (const auto& s: symbols_) {
                    if (SymEngine::eq(x, *s)) return true;
                }
                return false;
            }

            // Template method for `Symbol` like classes which do not have any
            // argument
            template<typename T> inline void remove_if_symbol_like(T& x)
            {
                result_ = condition_(x)
                    ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
            }

            // Template method for one argument function like classes
            template<typename Fun, typename Arg = const SymEngine::MatrixExpr, typename... Params>
            inline void remove_if_one_arg_f(
                Fun& x,
                std::function<SymEngine::RCP<Arg>()> get_arg,
                Params... params
            )
            {
                // We first check if the function will be removed
                if (condition_(x)) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                // Next we check if its argument will be removed
                else {
                    auto new_arg = apply(get_arg());
                    if (new_arg.is_null()) {
                        result_ = SymEngine::RCP<const SymEngine::Basic>();
                    }
                    else {
                        result_ = SymEngine::make_rcp<Fun>(
                            SymEngine::rcp_dynamic_cast<Arg>(new_arg), params...
                        );
                    }
                }
            }

        public:
            explicit RemoveVisitor(
                const SymEngine::vec_basic& symbols,
                std::function<bool(const SymEngine::Basic&)> condition
            ) : symbols_(symbols)
            {
                if (condition) {
                    condition_ = condition;
                }
                else {
                    condition_ = [=](const SymEngine::Basic& x) -> bool
                    {
                        return this->is_equal(x);
                    };
                }
            }

            inline SymEngine::RCP<const SymEngine::Basic> apply(
                const SymEngine::RCP<const SymEngine::Basic>& x
            )
            {
                if (condition_(*x)) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                } else {
                    x->accept(*this);
                }
                return result_;
            }

            void bvisit(const SymEngine::Basic& x);
            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const SymEngine::Integer& x);
            void bvisit(const SymEngine::Rational& x);
            void bvisit(const SymEngine::Complex& x);
            void bvisit(const SymEngine::Add& x);
            void bvisit(const SymEngine::Mul& x);
            void bvisit(const SymEngine::Constant& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::ZeroMatrix& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            void bvisit(const SymEngine::Trace& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixAdd& x);
            void bvisit(const SymEngine::MatrixMul& x);
            void bvisit(const SymEngine::MatrixDerivative& x);
    };

    // Remove given symbols from `x`
    inline SymEngine::RCP<const SymEngine::Basic> remove_if(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::vec_basic& symbols
    )
    {
        RemoveVisitor visitor(symbols);
        return visitor.apply(x);
    }
}
