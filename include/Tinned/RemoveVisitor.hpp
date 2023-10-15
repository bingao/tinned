/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of removal of specific symbols.

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

namespace Tinned
{
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

            void bvisit(const SymEngine::Basic& x);
            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
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
