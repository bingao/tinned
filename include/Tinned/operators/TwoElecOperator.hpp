/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of two-electron like operators.

   2023-10-27, Bin Gao:
   * rewrite member method get_args() that returns only the electronic state

   2023-09-22, Bin Gao:
   * first version
*/

#pragma once

#include <string>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_symbol.h>

#include "Tinned/perturbations/PertDependency.hpp"

#include "Tinned/operators/ElectronicState.hpp"

namespace Tinned
{
    // TwoElecOperator can be viewed as a tensor contraction of electron
    // repulsion integrals (ERI) and density matrix
    class TwoElecOperator: public SymEngine::MatrixSymbol
    {
        protected:
            // Electron state (may contain derivatives) that the operator
            // depends on
            SymEngine::RCP<const ElectronicState> state_;
            // dependencies_ stores perturbations that the operator depends on
            // and their maximum orders that can be differentiated
            PertDependency dependencies_;
            // derivatives_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivatives_;

        public:
            //! Constructor
            explicit TwoElecOperator(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const PertDependency& dependencies,
                const SymEngine::multiset_basic& derivatives = {}
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

    // Helper function to make a two-electron like operator
    inline SymEngine::RCP<const TwoElecOperator> make_2el_operator(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const PertDependency& dependencies = {}
    )
    {
        return SymEngine::make_rcp<const TwoElecOperator>(name, state, dependencies);
    }
}
