/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of keeping specific symbols while removing
   others.

   2023-10-19, Bin Gao:
   * moved from file Tinned/RemoveVisitor.hpp
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

#include "Tinned/RemoveVisitor.hpp"

namespace Tinned
{
    // Keeping symbols if they match any given ones while removing others
    class KeepVisitor: public SymEngine::BaseVisitor<KeepVisitor, RemoveVisitor>
    {
        protected:
            // Check inequality for `x` and symbols to be kept
            inline bool is_not_equal(const SymEngine::Basic& x) {
                for (const auto& s: symbols_) {
                    if (SymEngine::eq(x, *s)) return false;
                }
                return true;
            }

            // Template method for one argument function like classes
            template<typename Fun, typename Arg = const SymEngine::MatrixExpr, typename... Params>
            inline void remove_ifnot_one_arg_f(
                Fun& x,
                std::function<SymEngine::RCP<Arg>()> get_arg,
                Params... params
            )
            {
                // If the function will not be kept as whole, we then check if
                // its argument will be kept
                if (condition_(x)) {
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
                // The function will be kept as a whole
                else {
                    result_ = x.rcp_from_this();
                }
            }

        public:
            explicit KeepVisitor(
                const SymEngine::vec_basic& symbols
            ) : SymEngine::BaseVisitor<KeepVisitor, RemoveVisitor>(
                    symbols,
                    [=](const SymEngine::Basic& x) -> bool
                    {
                        return this->is_not_equal(x);
                    }
                )
            {
            }

            using RemoveVisitor::bvisit;
            //
            // Different from `RemoveVisitor`, the whole `Mul`, `MatrixMul`
            // and `HadamardProduct` will be kept whenever there is one factor
            // matches given symbols. Moreover, a function or an operator will
            // be kept if one of its argument matches given symbols.
            //
            void bvisit(const SymEngine::Mul& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            void bvisit(const SymEngine::Trace& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixMul& x);
    };

    // Keep given symbols in `x` while removing others
    inline SymEngine::RCP<const SymEngine::Basic> keep_if(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::vec_basic& symbols
    )
    {
        KeepVisitor visitor(symbols);
        return visitor.apply(x);
    }
}
