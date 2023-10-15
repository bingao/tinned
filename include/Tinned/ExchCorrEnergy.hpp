/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of exchange-correlation energy like
   functionals.

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
#include "Tinned/Perturbation.hpp"

namespace Tinned
{
    class ExchCorrEnergy: public SymEngine::FunctionWrapper
    {
        private:
            // derivative_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivative_;

        public:
            //! Constructor
            explicit ExchCorrEnergy(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const SymEngine::multiset_basic& derivative = {}
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

            // Get electronic state
            inline SymEngine::RCP<const ElectronicState> get_state() const
            {
                auto args = SymEngine::FunctionWrapper::get_args();
                return SymEngine::rcp_dynamic_cast<const ElectronicState>(args[0]);
            }

            // Get derivative
            inline SymEngine::multiset_basic get_derivative() const
            {
                return derivative_;
            }
    };
}
