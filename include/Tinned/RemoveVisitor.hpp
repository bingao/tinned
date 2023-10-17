/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of removal of specific symbols.

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
    // Rules for this removal procedure are that,
    // (1) Nothing left after removal, which is different from replacing
    //     symbols by zero.
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
    // (3) Symbols and their derivatives are different for the removal
    //     procedure.
    class RemoveVisitor: public SymEngine::BaseVisitor<RemoveVisitor>
    {
        protected:
            SymEngine::RCP<const SymEngine::Basic> result_;
            const SymEngine::vec_basic& symbols_;
            std::function<bool(const SymEngine::Basic&)> to_remove_;

            // Check equality for `x` and symbols to be removed
            inline bool eq_check(const SymEngine::Basic& x) {
                for (const auto& s: symbols_) {
                    if SymEngine::eq(x, *s) return true;
                }
                return false;
            }

            // Check inequality for `x` and symbols to be kept
            inline bool neq_check(const SymEngine::Basic& x) {
                for (const auto& s: symbols_) {
                    if SymEngine::eq(x, *s) return false;
                }
                return true;
            }

        public:
            explicit RemoveVisitor(
                const SymEngine::vec_basic& symbols,
                bool equivalence = true
            );

            inline SymEngine::RCP<const SymEngine::Basic> apply(
                const SymEngine::RCP<const SymEngine::Basic>& x
            )
            {
                if (to_remove(*x)) {
                    result_ = SymEngine::RCP();
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
    inline SymEngine::RCP<const SymEngine::Basic> eq_remove(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::vec_basic& symbols
    )
    {
        RemoveVisitor visitor(symbols, true);
        return visitor.apply(x);
    }

    // Remove symbols not in given ones from `x`
    inline SymEngine::RCP<const SymEngine::Basic> neq_remove(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::vec_basic& symbols
    )
    {
        RemoveVisitor visitor(symbols, false);
        return visitor.apply(x);
    }
}
