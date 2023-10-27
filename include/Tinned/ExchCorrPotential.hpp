/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of exchange-correlation (XC) potential like
   operators.

   2023-10-26, Bin Gao:
   * rewrite by using `CompositeFunction` and `DensityVector`, and built on top
     of `ExchCorrEnergy`

   2023-09-24, Bin Gao:
   * first version
*/

#pragma once

#include <string>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_symbol.h>

#include "Tinned/ElectronicState.hpp"
#include "Tinned/DensityVector.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/ExchCorrEnergy.hpp"

namespace Tinned
{
    class ExchCorrPotential: public SymEngine::MatrixSymbol
    {
        protected:
            SymEngine::RCP<const ElectronicState> state_;
            SymEngine::RCP<const NonElecFunction> weight_;
            SymEngine::RCP<const OneElecOperator> Omega_;

            // XC potential operator or its derivatives evaluated at grid
            // points
            SymEngine::RCP<const SymEngine::Basic> potential_;

        public:
            //! Constructor
            // `state`: electronic state like one-electron spin-orbital density matrix
            // `dependencies`: perturbation dependencies for the overlap distribution
            // `weight`: grid weight
            explicit ExchCorrPotential(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const PertDependency& dependencies,
                const SymEngine::RCP<const NonElecFunction>& weight
                    = SymEngine::make_rcp<const NonElecFunction>(
                        std::string("weight"),
                        PertDependency({})
                      )
            );
            // Constructor only for `diff_impl()`
            explicit ExchCorrPotential(
                const ExchCorrPotential& other,
                const SymEngine::RCP<const SymEngine::Symbol>& s
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get electronic state
            inline SymEngine::RCP<const ElectronicState> get_state() const
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
