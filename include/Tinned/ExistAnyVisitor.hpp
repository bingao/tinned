/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of checking if any given symbols exist in an
   expression.

   2024-10-08, Bin Gao:
   * first version
*/

#pragma once

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/mul.h>
#include <symengine/number.h>
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
    class ExistAnyVisitor: public SymEngine::BaseVisitor<ExistAnyVisitor>
    {
        protected:
            bool result_;
            SymEngine::set_basic symbols_;

            // Function template checking if `x` matches any given symbol
            template<typename T> inline bool exist_any_equivalence(T& x)
            {
                for (const auto& symbol: symbols_) {
                    if (symbol->__eq__(x)) {
                        result_ = true;
                        return true;
                    }
                }
                return false;
            }

            // Method called by objects to process their argument(s)
            inline bool apply_(const SymEngine::RCP<const SymEngine::Basic>& x)
            {
                x->accept(*this);
                return result_;
            }

        public:
            explicit ExistAnyVisitor(const SymEngine::set_basic& symbols)
                : symbols_(symbols) {}

            inline bool apply(const SymEngine::RCP<const SymEngine::Basic>& x)
            {
                result_ = false;
                return apply_(x);
            }

            void bvisit(const SymEngine::Basic& x);
            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const SymEngine::Number& x);
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

    // Helper function to check if any `symbols` exist in `x`
    inline bool exist_any(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::set_basic& symbols
    )
    {
        ExistAnyVisitor visitor(symbols);
        return visitor.apply(x);
    }
}
