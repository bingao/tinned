/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of finding a given symbol and all its
   differentiated ones.

   2023-10-28, Bin Gao:
   * first version
*/

#pragma once

#include <functional>
#include <set>

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
    class FindAllVisitor: public SymEngine::BaseVisitor<FindAllVisitor>
    {
        protected:
            SymEngine::set_basic result_;
            SymEngine::RCP<const SymEngine::Basic> symbol_;

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
            template<typename Fun, typename Arg>
            inline void remove_if_one_arg_f(
                Fun& x,
                const SymEngine::RCP<Arg>& arg,
                std::function<SymEngine::RCP<Fun>(const SymEngine::RCP<Arg>&)> constructor
            )
            {
                // We first check if the function will be removed
                if (condition_(x)) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                // Next we check if its argument will be removed
                else {
                    auto new_arg = apply(arg);
                    if (new_arg.is_null()) {
                        result_ = SymEngine::RCP<const SymEngine::Basic>();
                    }
                    else {
                        if (SymEngine::eq(*arg, *new_arg)) {
                            result_ = x.rcp_from_this();
                        }
                        else {
                            result_ = constructor(
                                SymEngine::rcp_dynamic_cast<Arg>(new_arg)
                            );
                        }
                    }
                }
            }

        public:
            explicit FindAllVisitor(
                const SymEngine::RCP<const SymEngine::Basic>& symbol,
                std::function<bool(const SymEngine::Basic&)> condition = {}
            ) : symbol_(symbol)
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

    // Find a given `symbol` and all its differentiated ones from `x`
    template<typename T>
    inline std::set<SymEngine::RCP<const T>, SymEngine::RCPBasicKeyLess> find_all(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::RCP<const SymEngine::Basic>& symbol
    )
    {
        std::set<SymEngine::RCP<const T>, SymEngine::RCPBasicKeyLess> result;
        FindAllVisitor visitor(symbol);
        for (auto& s: visitor.apply(x)) {
            if (SymEngine::is_a_sub<const T>(*s)) {
                result.insert(SymEngine::rcp_dynamic_cast<const T>(s));
            }
            else {
                throw SymEngine::SymEngineException(
                    "find_all() encounters an invalid type of: "+s->__str__()
                );
            }
        }
        return result;
    }
}
