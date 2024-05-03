/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of Lagrangian multipliers.

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

namespace Tinned
{
    class LagMultiplier: public SymEngine::MatrixSymbol
    {
        protected:
            // derivatives_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivatives_;

        public:
            //! Constructor
            // `derivatives` may be used only for `diff_impl()`
            explicit LagMultiplier(
                const std::string& name,
                const SymEngine::multiset_basic& derivatives = {}
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get derivatives
            inline SymEngine::multiset_basic get_derivatives() const
            {
                return derivatives_;
            }

            // Check if `x` is a same response parameter
            inline bool is_same_parameter(
                const SymEngine::RCP<const SymEngine::Basic>& x
            ) const
            {
                if (SymEngine::is_a_sub<const LagMultiplier>(*x)) {
                    auto op = SymEngine::rcp_dynamic_cast<const LagMultiplier>(x);
                    return get_name()==op->get_name() ? true : false;
                }
                return false;
            }
    };

    // Helper function to make a vector Lagrangian multipliers
    inline SymEngine::RCP<const LagMultiplier> make_lagrangian_multiplier(
        const std::string& name
    )
    {
        return SymEngine::make_rcp<const LagMultiplier>(name);
    }
}
