/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of exchange-correlation (XC) potential like
   operators.

   2024-05-04, Bin Gao:
   * add function `get_derivatives`

   2023-10-26, Bin Gao:
   * rewrite by using `CompositeFunction`

   2023-09-24, Bin Gao:
   * first version
*/

#pragma once

#include <map>
#include <string>
#include <set>
#include <utility>
#include <vector>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/symengine_assert.h>
#include <symengine/matrices/matrix_symbol.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/ExchCorrContraction.hpp"
#include "Tinned/FindAllVisitor.hpp"

namespace Tinned
{
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
            SymEngine::RCP<const SymEngine::MatrixExpr> potential_;

        public:
            //! Constructor
            // `state`: electronic state like one-electron spin-orbital density matrix
            // `Omega`: generalized overlap distribution vector
            // `weight`: grid weight
            explicit ExchCorrPotential(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const SymEngine::RCP<const OneElecOperator>& Omega,
                const SymEngine::RCP<const NonElecFunction>& weight
            );
            // Constructor only for `diff_impl()`
            explicit ExchCorrPotential(
                const ExchCorrPotential& other,
                const SymEngine::RCP<const SymEngine::Symbol>& s
            );
            // Constructor mainly used by different visitors
            explicit ExchCorrPotential(
                //const ExchCorrPotential& other,
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const SymEngine::RCP<const OneElecOperator>& Omega,
                const SymEngine::RCP<const NonElecFunction>& weight,
                const SymEngine::RCP<const SymEngine::MatrixExpr>& potential
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
            inline SymEngine::RCP<const SymEngine::MatrixExpr> get_potential() const
            {
                return potential_;
            }

            // Get all unique unperturbed and perturbed grid weights
            inline std::map<unsigned int, SymEngine::set_basic> get_weights() const
            {
                return find_all(potential_, get_weight());
            }

            // Get all unique unperturbed and perturbed electronic states
            inline std::map<unsigned int, SymEngine::set_basic> get_states() const
            {
                return find_all(potential_, get_state());
            }

            // Get derivatives, currently used only for `LaTeXifyVisitor` and
            // `OperatorEvaluator`
            inline SymEngine::multiset_basic get_derivatives() const
            {
                auto states = get_states();
                SYMENGINE_ASSERT(states.rbegin()->second.size()==1)
                auto state = states.rbegin()->second.begin();
                SYMENGINE_ASSERT(SymEngine::is_a_sub<const ElectronicState>(*(*state)))
                auto op = SymEngine::rcp_dynamic_cast<const ElectronicState>(*state);
                return op->get_derivatives();
            }

            // Get all unique unperturbed and perturbed generalized overlap
            // distribution vectors
            inline std::map<unsigned int, SymEngine::set_basic>
            get_overlap_distributions() const
            {
                return find_all(potential_, get_overlap_distribution());
            }

            // Get all unique orders of functional derivatives of XC energy density
            inline std::set<unsigned int> get_exc_orders() const
            {
                auto exc = find_all(potential_, make_exc_density(state_, Omega_, 1));
                std::set<unsigned int> orders;
                for (auto& e: exc) orders.insert(e.first);
                return orders;
            }

            // Get all terms in XC potential operator or its derivatives, each
            // is a product of (un)perturbed grid weight, XC functional
            // derivative vector, perturbed generalized density vectors and
            // (un)perturbed generalized overlap distribution.
            inline std::vector<SymEngine::RCP<const SymEngine::MatrixMul>>
            get_potential_terms() const
            {
                // Unperturbed case or when the generalized overlap
                // distribution does not depend on the applied perturbation(s)
                if (SymEngine::is_a_sub<const SymEngine::MatrixMul>(*potential_)) {
                    return std::vector<SymEngine::RCP<const SymEngine::MatrixMul>>({
                        SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(potential_)
                    });
                }
                // Perturbed case in particular the generalized overlap
                // distribution depends on the applied perturbation(s).
                // Constructor of `ExchCorrPotentail` has ensured the type of
                // `potential_` to be either `SymEngine::MatrixMul` or
                // `SymEngine::MatrixAdd`
                else {
                    std::vector<SymEngine::RCP<const SymEngine::MatrixMul>> terms;
                    auto potential = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixAdd>(potential_);
                    for (auto& contr: potential->get_args()) {
                        SYMENGINE_ASSERT(
                            SymEngine::is_a_sub<const SymEngine::MatrixMul>(*contr)
                        )
                        terms.push_back(
                            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(contr)
                        );
                    }
                    return terms;
                }
            }

            // Get all terms in XC potential operator or its derivatives, i.e.
            // (un)perturbed grid weights, XC functional derivative vectors,
            // perturbed generalized density vectors and (un)perturbed
            // generalized overlap distributions. Results are arranged in a
            // nested map. The key of the outermost map is (un)perturbed
            // generalized overlap distributions, whose value is
            // `ExcContractionMap`.
            inline VxcContractionMap get_potential_map() const
            {
                return extract_potential_map(potential_);
            }
    };

    // Helper function to make an XC potential like operator
    inline SymEngine::RCP<const ExchCorrPotential> make_xc_potential(
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
        return SymEngine::make_rcp<const ExchCorrPotential>(name, state, Omega, weight);
    }
}
