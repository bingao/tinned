/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of exchange-correlation potential like
   operators.

   2023-09-24, Bin Gao:
   * first version
*/

#pragma once

#include <string>

#include <symengine/dict.h>
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/matrix_symbol.h>

#include "Tinned/ElectronState.hpp"
#include "Tinned/Perturbation.hpp"

namespace Tinned
{
    class ExchCorrPotential: public SymEngine::MatrixSymbol
    {
        private:
            // Electron state that the operator depends on
            SymEngine::RCP<const ElectronState> state_;
            // derivative_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivative_;

        public:
            //! Constructor
            explicit ExchCorrPotential(
                const std::string& name,
                const SymEngine::RCP<const ElectronState>& state,
                const SymEngine::multiset_basic& derivative = {}
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::MatrixExpr> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get electron state
            inline SymEngine::RCP<const ElectronState> get_state() const
            {
                return state_;
            }

            // Get derivative
            inline SymEngine::multiset_basic get_derivative() const
            {
                return derivative_;
            }
    };
}
