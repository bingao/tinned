/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of state operators.

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

#include "Tinned/StateVector.hpp"

namespace Tinned
{
    class StateOperator: public SymEngine::MatrixSymbol
    {
        protected:
            // State vector
            SymEngine::RCP<const StateVector> state_;

        public:
            explicit StateOperator(
                const std::string& name,
                const SymEngine::RCP<const StateVector>& state
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get the state vector
            inline SymEngine::RCP<const StateVector> get_state() const
            {
                return state_;
            }

            // Get derivatives
            inline SymEngine::multiset_basic get_derivatives() const
            {
                return state_->get_derivatives();
            }
    };

    // Helper function to make a state operator
    inline SymEngine::RCP<const StateOperator> make_state_operator(
        const std::string& name,
        const SymEngine::RCP<const StateVector>& state
    )
    {
        return SymEngine::make_rcp<const StateOperator>(name, state);
    }
}