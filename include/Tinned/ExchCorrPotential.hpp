/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of exchange-correlation (XC) potential like
   operators.

   2023-10-26, Bin Gao:
   * rewrite by using `CompositeFunction` and `ExchCorrEnergy`

   2023-09-24, Bin Gao:
   * first version
*/

#pragma once

#include <string>
#include <set>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_symbol.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/ExchCorrEnergy.hpp"

namespace Tinned
{
    // Forward declaration for finding components in `potential_`
    class FindAllVisitor;

    // Make exchange-correlation potential, or the functional derivative of
    // exchange-correlation energy
    inline SymEngine::RCP<const ExchCorrEnergy> make_vxc(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const NonElecFunction>& weight
    )
    {
        return SymEngine::make_rcp<const ExchCorrEnergy>(name, state, Omega, weight, 1);
    }

    // Exchange-correlation (XC) potential like operators
    class ExchCorrPotential: public SymEngine::MatrixSymbol
    {
        protected:
            // Electronic state
            SymEngine::RCP<const ElectronicState> state_;
            // Overlap distribution
            SymEngine::RCP<const OneElecOperator> Omega_;
            // Grid weight
            SymEngine::RCP<const NonElecFunction> weight_;

            // XC potential operator or its derivatives evaluated at grid
            // points
            SymEngine::RCP<const SymEngine::Basic> potential_;

        public:
            //! Constructor
            // `state`: electronic state like one-electron spin-orbital density matrix
            // `Omega`: overlap distribution
            // `weight`: grid weight
            explicit ExchCorrPotential(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const SymEngine::RCP<const OneElecOperator>& Omega,
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
            // Constructor mainly used by different visitors
            explicit ExchCorrPotential(
                const ExchCorrPotential& other,
                const SymEngine::RCP<const SymEngine::Basic>& potential
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get grid weight
            inline SymEngine::RCP<const NonElecFunction> get_weight() const
            {
                return weight_;
            }

            // Get electronic state
            inline SymEngine::RCP<const ElectronicState> get_state() const
            {
                return state_;
            }

            // Get overlap distribution
            inline SymEngine::RCP<const OneElecOperator> get_overlap_distribution() const
            {
                return Omega_;
            }

            // Get XC potential operator or its derivatives evaluated at grid
            // points
            inline SymEngine::RCP<const SymEngine::Basic> get_potential() const
            {
                return potential_;
            }

            // Get all unique unperturbed and perturbed grid weights
            inline std::set<SymEngine::RCP<const NonElecFunction>,
                            SymEngine::RCPBasicKeyLess> get_weights() const
            {
                std::set<SymEngine::RCP<const NonElecFunction>,
                         SymEngine::RCPBasicKeyLess> weights;
                FindAllVisitor visitor(weight_);
                for (auto& w: visitor.apply(potential_)) {
                    weights.insert(
                        SymEngine::rcp_dynamic_cast<const NonElecFunction>(w)
                    );
                }
                return weights;
            }

            // Get all unique unperturbed and perturbed electronic states
            inline std::set<SymEngine::RCP<const ElectronicState>,
                            SymEngine::RCPBasicKeyLess> get_states() const
            {
                std::set<SymEngine::RCP<const ElectronicState>,
                         SymEngine::RCPBasicKeyLess> states;
                FindAllVisitor visitor(state_);
                for (auto& s: visitor.apply(potential_)) {
                    states.insert(
                        SymEngine::rcp_dynamic_cast<const ElectronicState>(s)
                    );
                }
                return states;
            }

            // Get all unique unperturbed and perturbed overlap distributions
            inline std::set<SymEngine::RCP<const OneElecOperator>,
                            SymEngine::RCPBasicKeyLess> get_overlap_distributions() const
            {
                std::set<SymEngine::RCP<const OneElecOperator>,
                         SymEngine::RCPBasicKeyLess> Omegas;
                FindAllVisitor visitor(Omega_);
                for (auto& w: visitor.apply(potential_)) {
                    Omegas.insert(
                        SymEngine::rcp_dynamic_cast<const OneElecOperator>(w)
                    );
                }
                return Omegas;
            }

            // Get all unique orders of functional derivatives of XC energy density
            inline std::set<unsigned int> get_exc_orders() const
            {
                std::set<unsigned int> orders;
                FindAllVisitor visitor(make_exc_density(state_, Omega_, 1));
                for (auto& exc: visitor.apply(potential_)) {
                    orders.insert(
                        SymEngine::rcp_dynamic_cast<const CompositeFunction>(exc)->get_order()
                    );
                }
                return orders;
            }
    };
}
