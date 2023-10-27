/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of exchange-correlation (XC) energy like
   functionals.

   2023-10-25, Bin Gao:
   * rewrite by using `CompositeFunction`

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
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/NonElecFunction.hpp"

namespace Tinned
{
    // Forward declaration for finding components in `energy_`
    class FindAllVisitor;

    // Make generalized density vector
    //FIXME: change `ElectronicState` to OneElecDensity?
    inline SymEngine::RCP<const SymEngine::Basic> make_density_vector(
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega
    )
    {
        return SymEngine::trace(SymEngine::matrix_mul({Omega, state}));
    }

    // Make XC energy density
    inline SymEngine::RCP<const CompositeFunction> make_exc_density(
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const unsigned int order
    )
    {
        return SymEngine::make_rcp<const CompositeFunction>(
            // XC energy density as the outer function
            std::string("exc"),
            // Generalized density vector as the inner function
            make_density_vector(state, Omega),
            order
        );
    }

    // Exchange-correlation (XC) energy like functionals
    class ExchCorrEnergy: public SymEngine::FunctionWrapper
    {
        protected:
            // XC energy or its derivatives evaluated at grid points
            SymEngine::RCP<const SymEngine::Basic> energy_;

        public:
            //! Constructor
            // `state`: electronic state like one-electron spin-orbital density matrix
            // `Omega`: overlap distribution
            // `weight`: grid weight
            // `order`: order of functional derivatives of XC energy density
            explicit ExchCorrEnergy(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const SymEngine::RCP<const OneElecOperator>& Omega,
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

            // Get grid weight
            inline SymEngine::RCP<const SymEngine::Basic> get_weight() const
            {
                return get_args()[0];
            }

            // Get electronic state
            inline SymEngine::RCP<const SymEngine::Basic> get_state() const
            {
                return get_args()[1];
            }

            // Get overlap distribution
            inline SymEngine::RCP<const SymEngine::Basic> get_overlap_distribution() const
            {
                return get_args()[2];
            }

            // Get XC energy or its derivatives evaluated at grid points
            inline SymEngine::RCP<const SymEngine::Basic> get_energy() const
            {
                return energy_;
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

            // Get all unique unperturbed and perturbed overlap distributions
            inline std::set<SymEngine::RCP<const OneElecOperator>,
                            SymEngine::RCPBasicKeyLess> get_overlap_distributions() const
            {
                std::set<SymEngine::RCP<const OneElecOperator>,
                         SymEngine::RCPBasicKeyLess> Omegas;
                FindAllVisitor visitor(get_overlap_distribution());
                for (auto& w: visitor.apply(energy_)) {
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
                FindAllVisitor visitor(
                    make_exc_density(get_state(), get_overlap_distribution(), 0)
                );
                for (auto& exc: visitor.apply(energy_)) {
                    orders.insert(
                        SymEngine::rcp_dynamic_cast<const CompositeFunction>(exc)->get_order()
                    );
                }
                return orders;
            }
    };
}
