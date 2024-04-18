/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of exponential maps.

   2024-04-17, Bin Gao:
   * first version
*/

#pragma once

#include <string>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_symbol.h>

#include "Tinned/StateOperator.hpp"

namespace Tinned
{
    class ExpAdjointHamiltonian: public SymEngine::MatrixSymbol
    {
        protected:
            // State operator
            SymEngine::RCP<const StateOperator> state_operator_;
            // Hamiltonian
            SymEngine::RCP<const SymEngine::Basic> hamiltonian_;

        public:
            explicit ExpAdjointHamiltonian(
                const std::string& name,
                const SymEngine::RCP<const StateOperator>& state_operator,
                const SymEngine::RCP<const SymEngine::Basic>& hamiltonian
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get the state operator
            inline SymEngine::RCP<const StateOperator> get_state_operator() const
            {
                return state_operator_;
            }

            // Get the Hamiltonian
            inline SymEngine::RCP<const SymEngine::Basic> get_hamiltonian() const
            {
                return hamiltonian_;
            }
    };

    // Helper function to make an exponential map
    inline SymEngine::RCP<const ExpAdjointHamiltonian> make_eadj_hamiltonian(
        const std::string& name,
        const SymEngine::RCP<const StateOperator>& state_operator,
        const SymEngine::RCP<const SymEngine::Basic>& hamiltonian
    )
    {
        return SymEngine::make_rcp<const ExpAdjointHamiltonian>(
            name, state_operator, hamiltonian
        );
    }
}
