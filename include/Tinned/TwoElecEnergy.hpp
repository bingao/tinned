/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of two-electron like energies.

   2023-11-22, Bin Gao:
   * first version
*/

#pragma once

#include <string>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/number.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

#include "Tinned/ElectronicState.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/TwoElecOperator.hpp"

namespace Tinned
{
    class TwoElecEnergy: public SymEngine::FunctionWrapper
    {
        protected:
            // Electron states (may contain derivatives) that the two-electron
            // energy is represented as Tr[G(`inner_`)*`outer_`]
            SymEngine::RCP<const ElectronicState> inner_;
            SymEngine::RCP<const ElectronicState> outer_;
            // dependencies_ stores perturbations that the operator depends on
            // and their maximum orders that can be differentiated
            PertDependency dependencies_;
            // derivatives_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivatives_;

        public:
            //! Constructor
            // `derivatives` may only be used for `diff_impl()`
            explicit TwoElecEnergy(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& inner,
                const SymEngine::RCP<const ElectronicState>& outer,
                const PertDependency& dependencies,
                const SymEngine::multiset_basic& derivatives = {}
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            //SymEngine::vec_basic get_args() const override;

            SymEngine::RCP<const SymEngine::Basic> create(
                const SymEngine::vec_basic &v
            ) const override;
            SymEngine::RCP<const SymEngine::Number> eval(long bits) const override;
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol> &s
            ) const override;

            // Get inner electronic state
            inline SymEngine::RCP<const ElectronicState> get_inner_state() const
            {
                return inner_;
            }

            // Get outer electronic state
            inline SymEngine::RCP<const ElectronicState> get_outer_state() const
            {
                return outer_;
            }

            // Get dependencies
            inline PertDependency get_dependencies() const
            {
                return dependencies_;
            }

            // Get derivatives
            inline SymEngine::multiset_basic get_derivatives() const
            {
                return derivatives_;
            }
    };

    // Helper function to make two-electron like energies from
    // `ElectronicState`
    inline SymEngine::RCP<const TwoElecEnergy> make_2el_energy(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& inner,
        const SymEngine::RCP<const ElectronicState>& outer,
        const PertDependency& dependencies = {}
    )
    {
        return SymEngine::make_rcp<const TwoElecEnergy>(
            name, inner, outer, dependencies
        );
    }

    // Helper function to make two-electron like energies from
    // `TwoElecOperator`
    inline SymEngine::RCP<const TwoElecEnergy> make_2el_energy(
        const std::string& name,
        const SymEngine::RCP<const TwoElecOperator>& G
    )
    {
        return SymEngine::make_rcp<const TwoElecEnergy>(
            name, G->get_state(), G->get_state(), G->get_dependencies()
        );
    }
}
