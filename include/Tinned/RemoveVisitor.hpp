/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of removal of specific symbols.

   2023-10-18, Bin Gao:
   * two classes RemoveIfVisitor and RemoveIfNotVisitor for removing matching
     symbols and keeping matching symbols while removing others

   2023-10-16, Bin Gao:
   * formal rules for the removal procedure

   2023-09-25, Bin Gao:
   * first version
*/

#pragma once

#include <functional>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

namespace Tinned
{
    // Removing symbols if they match any given ones
    //
    // Rules for the removal procedure are that,
    //
    // (1) Nothing left after removal, which is different from replacing
    //     symbols by zero.
    //
    // (2) If a symbol is involved in an operation, the operation should still
    //     be valid/meaning after removal. Otherwise, either the whole
    //     operation will be removed or such a removal is not allowed.
    //
    //     For example, a symbol removed from an addition or a multiplication
    //     is allowed and is the same as replacing it by zero. That is, the
    //     multiplication will be removed if one of its factor is removed.
    //
    //     Removal of an exponent from a power, removal of a variable that a
    //     derivative is with respect to, and removal of an element from a
    //     matrix, all are not allowed, because the left power and
    //     differentiation operations, and the matrix become invalid.
    //
    // (3) Symbols and their derivatives are different for the removal
    //     procedure.
    class RemoveIfVisitor: public SymEngine::BaseVisitor<RemoveIfVisitor>
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
            explicit RemoveIfVisitor(
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

    // Removing symbols if they do not match any given ones, or more exactly,
    // keeping symbols if they match any given ones
    class RemoveIfNotVisitor: public SymEngine::BaseVisitor<RemoveIfNotVisitor, RemoveIfVisitor>
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
            explicit RemoveIfNotVisitor(
                const SymEngine::vec_basic& symbols
            ) : SymEngine::BaseVisitor<RemoveIfNotVisitor, RemoveIfVisitor>(
                    symbols,
                    [=](const SymEngine::Basic& x) -> bool
                    {
                        return this->is_not_equal(x);
                    }
                )
            {
            }

            using RemoveIfVisitor::bvisit;
            //
            // Different from `RemoveIfVisitor`, the whole `Mul`, `MatrixMul`
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

    // Remove given symbols from `x`
    inline SymEngine::RCP<const SymEngine::Basic> remove_if(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::vec_basic& symbols
    )
    {
        RemoveIfVisitor visitor(symbols);
        return visitor.apply(x);
    }

    // Remove symbols if not in given ones from `x`, or more exactly, keep
    // given symbols in `x`
    inline SymEngine::RCP<const SymEngine::Basic> remove_if_not(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::vec_basic& symbols
    )
    {
        RemoveIfNotVisitor visitor(symbols);
        return visitor.apply(x);
    }
}
