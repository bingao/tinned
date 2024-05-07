/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of finding a given symbol and all its
   differentiated ones.

   2024-05-07, Bin Gao:
   * change `find_all` to a function that returns `SymEngine::set_basic`

   2023-10-28, Bin Gao:
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

#include "Tinned/PertDependency.hpp"

namespace Tinned
{
    class FindAllVisitor: public SymEngine::BaseVisitor<FindAllVisitor>
    {
        protected:
            SymEngine::set_basic result_;
            SymEngine::RCP<const SymEngine::Basic> symbol_;

            // Template method for objects that compares only the names
            template<typename T> inline void find_only_name(T& x)
            {
                if (SymEngine::is_a_sub<T>(*symbol_)) {
                    auto s = SymEngine::rcp_dynamic_cast<T>(symbol_);
                    if (x.get_name()==s->get_name()) result_.insert(x.rcp_from_this());
                }
            }

            // Template method for objects that requires equivalence comparison
            template<typename T> inline bool find_equivalence(T& x)
            {
                if (symbol_->__eq__(x)) {
                    result_.insert(x.rcp_from_this());
                    return true;
                }
                return false;
            }

            // Template method for objects without arguments, we compare only
            // their names and dependencies
            template<typename T> inline bool find_with_dependencies(T& x)
            {
                if (SymEngine::is_a_sub<T>(*symbol_)) {
                    auto s = SymEngine::rcp_dynamic_cast<T>(symbol_);
                    if (x.get_name()==s->get_name() &&
                        eq_dependency(x.get_dependencies(), s->get_dependencies())) {
                        result_.insert(x.rcp_from_this());
                        return true;
                    }
                }
                return false;
            }

        public:
            explicit FindAllVisitor(
                const SymEngine::RCP<const SymEngine::Basic>& symbol
            ) : symbol_(symbol) {}

            inline SymEngine::set_basic apply(
                const SymEngine::RCP<const SymEngine::Basic>& x
            )
            {
                x->accept(*this);
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

    // Helper function to find a given `symbol` and all its differentiated ones in `x`
    inline SymEngine::set_basic find_all(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::RCP<const SymEngine::Basic>& symbol
    )
    {
        FindAllVisitor visitor(symbol);
        return visitor.apply(x);
    }
}
