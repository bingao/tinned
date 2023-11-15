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
#include <set>
#include <vector>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/number.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/symengine_rcp.h>
#include <symengine/symengine_assert.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/ExchCorrContraction.hpp"
#include "Tinned/FindAllVisitor.hpp"

namespace Tinned
{
    // Make generalized density vector
    //FIXME: change `ElectronicState` to OneElecDensity?
    inline SymEngine::RCP<const SymEngine::Basic> make_density_vector(
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega
    )
    {
        return SymEngine::trace(SymEngine::matrix_mul({Omega, state}));
    }

    // Make XC energy density from electronic state, overlap distribution and
    // the order of XC energy functional
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
            // `Omega`: generalized overlap distribution vector
            // `weight`: grid weight
            // `order`: order of functional derivatives of XC energy density
            explicit ExchCorrEnergy(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const SymEngine::RCP<const OneElecOperator>& Omega,
                const SymEngine::RCP<const NonElecFunction>& weight,
                const unsigned int order = 0
            );
            // Constructor only for `diff_impl()`
            explicit ExchCorrEnergy(
                const ExchCorrEnergy& other,
                const SymEngine::RCP<const SymEngine::Symbol>& s
            );
            // Constructor mainly used by different visitors
            explicit ExchCorrEnergy(
                const ExchCorrEnergy& other,
                const SymEngine::RCP<const SymEngine::Basic>& energy
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
            inline SameTypeSet<const NonElecFunction> get_weights() const
            {
                return find_all<NonElecFunction>(energy_, get_weight());
            }

            // Get all unique unperturbed and perturbed electronic states
            inline SameTypeSet<const ElectronicState> get_states() const
            {
                return find_all<ElectronicState>(energy_, get_state());
            }

            // Get all unique unperturbed and perturbed generalized overlap
            // distribution vectors
            inline SameTypeSet<const OneElecOperator> get_overlap_distributions() const
            {
                return find_all<OneElecOperator>(energy_, get_overlap_distribution());
            }

            // Get all unique orders of functional derivatives of XC energy density
            inline std::set<unsigned int> get_exc_orders() const
            {
                auto exc = find_all<CompositeFunction>(
                    energy_,
                    make_exc_density(
                        SymEngine::rcp_dynamic_cast<const ElectronicState>(
                            get_state()
                        ),
                        SymEngine::rcp_dynamic_cast<const OneElecOperator>(
                            get_overlap_distribution()
                        ),
                        0
                    )
                );
                std::set<unsigned int> orders;
                for (auto& e: exc) orders.insert(e->get_order());
                return orders;
            }

            // Get all terms in XC energy or its derivatives, each is a product
            // of (un)perturbed weight, XC functional derivative vector and
            // perturbed generalized density vectors
            inline std::vector<SymEngine::RCP<const SymEngine::Mul>>
            get_energy_terms() const
            {
                // Unperturbed or the first-order case
                if (SymEngine::is_a_sub<const SymEngine::Mul>(*energy_)) {
                    return std::vector<SymEngine::RCP<const SymEngine::Mul>>({
                        SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(energy_)
                    });
                }
                // Perturbed case, constructor of `ExchCorrEnergy` has ensured
                // the type of `energy_` to be either `SymEngine::Mul` or
                // `SymEngine::Add`
                else {
                    std::vector<SymEngine::RCP<const SymEngine::Mul>> terms;
                    auto energy = SymEngine::rcp_dynamic_cast<const SymEngine::Add>(energy_);
                    auto contractions = energy->get_args();
                    // No coefficient exists in the XC energy derivatives
                    SYMENGINE_ASSERT(
                        !SymEngine::is_a_sub<const SymEngine::Number>(*contractions.front())
                    )
                    for (const auto& contr: contractions) {
                        SYMENGINE_ASSERT(
                            SymEngine::is_a_sub<const SymEngine::Mul>(*contr)
                        )
                        terms.push_back(
                            SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(contr)
                        );
                    }
                    return terms;
                }
            }

            // Get all terms in XC energy or its derivatives, i.e.
            // (un)perturbed weights, XC functional derivative vectors and
            // perturbed generalized density vectors. Results are arranged in a
            // nested map. The key of the outer map is (un)perturbed weights,
            // and the value of the outer map is another map whose key is the
            // order of functional derivatives of XC energy density, and whose
            // value is the corresponding perturbed generalized density
            // vectors.
            inline ExcContractionMap get_energy_map() const
            {
                // Unperturbed or the first-order case
                if (SymEngine::is_a_sub<const SymEngine::Mul>(*energy_)) {
                    auto contr_term = extract_exc_contraction(
                        SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(energy_)
                    );
                    return ExcContractionMap({
                        {
                            std::get<0>(contr_term),
                            ExcDensityContractionMap({
                                {std::get<1>(contr_term), std::get<2>(contr_term)}
                            })
                        }
                    });
                }
                // Perturbed case, constructor of `ExchCorrEnergy` has ensured
                // the type of `energy_` to be either `SymEngine::Mul` or
                // `SymEngine::Add`
                else {
                    ExcContractionMap contr_map;
                    auto energy = SymEngine::rcp_dynamic_cast<const SymEngine::Add>(energy_);
                    auto contractions = energy->get_args();
                    // No coefficient exists in the XC energy derivatives
                    SYMENGINE_ASSERT(
                        !SymEngine::is_a_sub<const SymEngine::Number>(*contractions.front())
                    )
                    for (const auto& contr: contractions) {
                        SYMENGINE_ASSERT(
                            SymEngine::is_a_sub<const SymEngine::Mul>(*contr)
                        )
                        auto contr_term = extract_exc_contraction(
                            SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(contr)
                        );
                        add_exc_contraction(contr_map, contr_term);
                    }
                    return contr_map;
                }
            }
    };

    // Helper function to make XC energy like functionals
    inline SymEngine::RCP<const ExchCorrEnergy> make_xc_energy(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const NonElecFunction>& weight
            = SymEngine::make_rcp<const NonElecFunction>(
                std::string("weight"),
                PertDependency({})
              )
    )
    {
        return SymEngine::make_rcp<const ExchCorrEnergy>(name, state, Omega, weight);
    }
}
