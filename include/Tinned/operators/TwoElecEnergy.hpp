/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of two-electron like energies.

   2024-05-11, Bin Gao:
   * store SymEngine::RCP<const TwoElecOperator> `G_(inner_)` and
     SymEngine::RCP<const ElectronicState> `outer_` as member variables

   2024-05-09, Bin Gao:
   * sort states according to their `compare` function

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

#include "Tinned/perturbations/PertDependency.hpp"

#include "Tinned/operators/ElectronicState.hpp"
#include "Tinned/operators/TwoElecOperator.hpp"

namespace Tinned
{
    class TwoElecEnergy: public SymEngine::FunctionWrapper
    {
        protected:
            // Two-electron energy is represented as tr[\frac{1}{2}*`G_(inner_)`*`outer_`]
            SymEngine::RCP<const TwoElecOperator> G_;
            SymEngine::RCP<const ElectronicState> outer_;

        public:
            //! Constructor for given density matrices, the two-electron
            //  operator will be constructed from `name`, `inner`,
            //  `dependencies` and `derivatives`
            explicit TwoElecEnergy(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& inner,
                const SymEngine::RCP<const ElectronicState>& outer,
                const PertDependency& dependencies,
                const SymEngine::multiset_basic& derivatives = {}
            );

            //! Constructor for a given two-electron operator and an outer
            //  density matrix
            explicit TwoElecEnergy(
                const SymEngine::RCP<const TwoElecOperator>& G,
                const SymEngine::RCP<const ElectronicState>& outer
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            SymEngine::RCP<const SymEngine::Basic> create(
                const SymEngine::vec_basic &v
            ) const override;
            SymEngine::RCP<const SymEngine::Number> eval(long bits) const override;
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol> &s
            ) const override;

            // Get two-electron operator
            inline SymEngine::RCP<const TwoElecOperator> get_2el_operator() const
            {
                return G_;
            }

            // Get inner electronic state
            inline SymEngine::RCP<const ElectronicState> get_inner_state() const
            {
                return G_->get_state();
            }

            // Get outer electronic state
            inline SymEngine::RCP<const ElectronicState> get_outer_state() const
            {
                return outer_;
            }

            // Get dependencies
            inline PertDependency get_dependencies() const
            {
                return G_->get_dependencies();
            }

            // Get derivatives
            inline SymEngine::multiset_basic get_derivatives() const
            {
                return G_->get_derivatives();
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

    // Helper function to make a two-electron like energy from
    // `TwoElecOperator`
    inline SymEngine::RCP<const TwoElecEnergy> make_2el_energy(
        const SymEngine::RCP<const TwoElecOperator>& G,
        const SymEngine::RCP<const ElectronicState>& outer = SymEngine::RCP<const ElectronicState>()
    )
    {
        return outer.is_null()
            ? SymEngine::make_rcp<const TwoElecEnergy>(G, G->get_state())
            : SymEngine::make_rcp<const TwoElecEnergy>(G, outer);
    }
}
