/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of exchange-correlation (XC) energy like
   functionals.

   2023-10-25, Bin Gao:
   * rewrite by using `CompositeFunction` and `DensityVector`

   2023-09-23, Bin Gao:
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
#include "Tinned/DensityVector.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/PertDependency.hpp"

namespace Tinned
{
    // Forward declaration for finding components in `energy_`
    class FindAllVisitor;

    class ExchCorrEnergy: public SymEngine::FunctionWrapper
    {
        protected:
            // Get XC energy density
            inline SymEngine::RCP<const CompositeFunction> get_exc_density(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const unsigned int order
            ) const
            {
                return SymEngine::make_rcp<const CompositeFunction>(
                    // XC energy density as the outer function
                    name,
                    // Generalized density vector as the inner function
                    SymEngine::make_rcp<const DensityVector>(
                        std::string("rho"),
                        state
                    ),
                    order
                );
            }

            // XC energy or its derivatives evaluated at grid points
            SymEngine::RCP<const SymEngine::Basic> energy_;

        public:
            //! Constructor
            // `state`: electronic state like one-electron spin-orbital density matrix
            // `weight`: grid weight
            // `order`: order of functional derivatives of XC energy density
            explicit ExchCorrEnergy(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const SymEngine::RCP<const NonElecFunction>& weight
                    = SymEngine::make_rcp<const NonElecFunction>(
                        std::string("weight"),
                        PertDependency({})
                      ),
                const unsigned int order = 0
            );
            // Constructor only for `diff_impl()`
            explicit ExchCorrEnergy(
                const ExchCorrEnergy& other,
                const SymEngine::RCP<const SymEngine::Symbol>& s
                //const SymEngine::RCP<const SymEngine::Basic> energy
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;

            SymEngine::RCP<const SymEngine::Basic> create(
                const SymEngine::vec_basic &v
            ) const override;
            SymEngine::RCP<const SymEngine::Number> eval(long bits) const override;
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol> &s
            ) const override;

            // Get electronic state
            inline SymEngine::RCP<const SymEngine::Basic> get_state() const
            {
                return get_args()[0];
            }

            // Get grid weight
            inline SymEngine::RCP<const SymEngine::Basic> get_weight() const
            {
                return get_args()[1];
            }

            // Get XC energy or its derivatives evaluated at grid points
            inline SymEngine::RCP<const SymEngine::Basic> get_energy() const
            {
                return energy_;
            }

            // Get all unique unperturbed and perturbed electronic states
            inline std::set<SymEngine::RCP<const ElectronicState>,
                            SymEngine::RCPBasicKeyLess> get_states() const
            {
                std::set<SymEngine::RCP<const ElectronicState>,
                         SymEngine::RCPBasicKeyLess> states;
                FindAllVisitor visitor(get_state());
                for (auto& s: visitor.apply(energy_)) {
                    states.insert(
                        SymEngine::rcp_dynamic_cast<const ElectronicState>(s)
                    );
                }
                return states;
            }

            // Get all unique unperturbed and perturbed grid weights
            inline std::set<SymEngine::RCP<const NonElecFunction>,
                            SymEngine::RCPBasicKeyLess> get_weights() const
            {
                std::set<SymEngine::RCP<const NonElecFunction>,
                         SymEngine::RCPBasicKeyLess> weights;
                FindAllVisitor visitor(get_weight());
                for (auto& w: visitor.apply(energy_)) {
                    weights.insert(
                        SymEngine::rcp_dynamic_cast<const NonElecFunction>(w)
                    );
                }
                return weights;
            }

            // Get all unique orders of functional derivatives of XC energy density
            inline std::set<unsigned int> get_exc_orders() const
            {
                std::set<unsigned int> orders;
                FindAllVisitor visitor(get_exc_density(get_name(), get_state(), 0));
                for (auto& exc: visitor.apply(energy_)) {
                    orders.insert(
                        SymEngine::rcp_dynamic_cast<const CompositeFunction>(exc)->get_order()
                    );
                }
                return orders;
            }
    };
}
