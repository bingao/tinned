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

#include <map>
#include <string>
#include <set>
#include <tuple>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/number.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_exception.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/KeepVisitor.hpp"
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
            // Extract grid weight, XC functional derivative and generalized
            // density vectors from their contraction
            inline std::tuple<SymEngine::RCP<const NonElecFunction>,
                              unsigned int,
                              SymEngine::RCP<const SymEngine::Basic>>
            extract_exc_contraction(
                const SymEngine::RCP<const SymEngine::Mul>& contraction
            ) const
            {
                SymEngine::RCP<const NonElecFunction> weight;
                SymEngine::RCP<const CompositeFunction> exc;
                unsigned int order;
                SymEngine::vec_basic factors;
                bool found_weight = false;
                bool found_order = false;
                bool found_factors = false;
                for (const auto& arg: contraction->get_args()) {
                    // Coefficient
                    if (SymEngine::is_a_sub<const SymEngine::Number>(*arg)) {
                        factors.push_back(arg);
                    }
                    // Grid weight
                    else if (SymEngine::is_a_sub<const NonElecFunction>(*arg)) {
                        if (found_weight) throw SymEngine::SymEngineException(
                            "Two grid weights got from the XC energy density contraction "
                            + contraction->__str__()
                        );
                        weight = SymEngine::rcp_dynamic_cast<const NonElecFunction>(arg);
                        found_weight = true;
                    }
                    // XC functional derivative
                    else if (SymEngine::is_a_sub<const CompositeFunction>(*arg)) {
                        if (found_order) throw SymEngine::SymEngineException(
                            "Two XC functional derivatives got from the XC energy density contraction "
                            + contraction->__str__()
                        );
                        exc = SymEngine::rcp_dynamic_cast<const CompositeFunction>(arg);
                        order = exc->get_order();
                        found_order = true;
                    }
                    // Generalized density vector, sum of generalized density
                    // vectors, or power of (sum of) generalized density
                    // vector(s)
                    else if (
                        SymEngine::is_a_sub<const SymEngine::Trace>(*arg) ||
                        SymEngine::is_a_sub<const SymEngine::Add>(*arg) ||
                        SymEngine::is_a_sub<const SymEngine::Pow>(*arg)
                    ) {
                        factors.push_back(arg);
                        found_factors = true;
                    }
                    else {
                        throw SymEngine::SymEngineException(
                            "Invalid type from the XC energy density contraction "
                            + contraction->__str__()
                        );
                    }
                }
                if (found_weight && found_order && (found_factors || order == 0)) {
                    // For unperturbed case, there is no contractions of the
                    // functional derivative vectors with the perturbed
                    // generalized density vectors
                    if (order == 0) {
                        if (found_factors) throw SymEngine::SymEngineException(
                            "Invalid factors for unperturbed XC energy functional "
                            + contraction->__str__()
                        );
                        return std::make_tuple(
                            weight,
                            0,
                            // Users should be aware that this is a null
                            // density vector by using `is_null()` for
                            // unperturbed case (`order` is 0)
                            SymEngine::RCP<const SymEngine::Basic>()
                        );
                    }
                    else {
                        auto dens_vectors = SymEngine::mul(factors);
                        // Last, we need to make sure there is neither grid
                        // weights nor XC functional derivatives in the
                        // generalized density vectors we got
                        auto factors_left = keep_if(
                            dens_vectors, SymEngine::set_basic({weight, exc})
                        );
                        if (factors_left.is_null()) {
                            return std::make_tuple(weight, order, dens_vectors);
                        }
                        else {
                            throw SymEngine::SymEngineException(
                                "Invalid grid weights and/or XC functional derivatives "
                                + factors_left->__str__()
                                + " in "
                                + dens_vectors->__str__()
                            );
                        }
                    }
                }
                else {
                    throw SymEngine::SymEngineException(
                        "Terms missing ("
                        + std::to_string(found_weight) + "/"
                        + std::to_string(found_order) + "/"
                        + std::to_string(found_factors)
                        + ") from the XC energy density contraction "
                        + contraction->__str__()
                    );
                }
            }

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
            inline std::set<SymEngine::RCP<const NonElecFunction>,
                            SymEngine::RCPBasicKeyLess> get_weights() const
            {
                return find_all<NonElecFunction>(energy_, get_weight());
            }

            // Get all unique unperturbed and perturbed electronic states
            inline std::set<SymEngine::RCP<const ElectronicState>,
                            SymEngine::RCPBasicKeyLess> get_states() const
            {
                return find_all<ElectronicState>(energy_, get_state());
            }

            // Get all unique unperturbed and perturbed overlap distributions
            inline std::set<SymEngine::RCP<const OneElecOperator>,
                            SymEngine::RCPBasicKeyLess> get_overlap_distributions() const
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

            // Get all terms in XC energy or its derivatives, i.e.
            // (un)perturbed weights, perturbed generalized density vectors and
            // the XC functional derivative vectors. Results are arranged in a
            // nested map. The key of the outer map is (un)perturbed weights,
            // and the value of the outer map is still a map whose key is the
            // order of functional derivatives of XC energy density, and value
            // is the corresponding perturbed generalized density vectors.
            inline std::map<SymEngine::RCP<const NonElecFunction>,
                            std::map<unsigned int, SymEngine::RCP<const SymEngine::Basic>>,
                            SymEngine::RCPBasicKeyLess>
            get_energy_terms() const
            {
                // Unperturbed or first-order case
                if (SymEngine::is_a_sub<const SymEngine::Mul>(*energy_)) {
                    auto terms = extract_exc_contraction(
                        SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(energy_)
                    );
                    return std::map<SymEngine::RCP<const NonElecFunction>,
                                    std::map<unsigned int, SymEngine::RCP<const SymEngine::Basic>>,
                                    SymEngine::RCPBasicKeyLess>({
                        {
                            std::get<0>(terms),
                            std::map<unsigned int, SymEngine::RCP<const SymEngine::Basic>>({
                                {std::get<1>(terms), std::get<2>(terms)}
                            })
                        }
                    });
                }
                // Perturbed case, constructor of `ExchCorrEnergy` has ensured
                // the type of `energy_` to be either `SymEngine::Mul` or
                // `SymEngine::Add`
                else {
                    std::map<SymEngine::RCP<const NonElecFunction>,
                             std::map<unsigned int, SymEngine::RCP<const SymEngine::Basic>>,
                             SymEngine::RCPBasicKeyLess> terms;
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
                        auto contr_terms = extract_exc_contraction(
                            SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(contr)
                        );
                        // Check if the grid weight exists in the outer map
                        auto iter = terms.find(std::get<0>(contr_terms));
                        if (iter == terms.end()) {
                            terms.emplace(
                                std::get<0>(contr_terms),
                                std::map<unsigned int, SymEngine::RCP<const SymEngine::Basic>>({
                                    {std::get<1>(contr_terms), std::get<2>(contr_terms)}
                                })
                            );
                        }
                        else {
                            // Check if the functional derivative of XC energy
                            // density exists in the inner map
                            auto jter = iter->second.find(std::get<1>(contr_terms));
                            if (jter == iter->second.end()) {
                                iter->second.emplace(
                                    std::get<1>(contr_terms), std::get<2>(contr_terms)
                                );
                            }
                            else {
                                jter->second = SymEngine::add(
                                    jter->second, std::get<2>(contr_terms)
                                );
                            }
                        }
                    }
                    return terms;
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
