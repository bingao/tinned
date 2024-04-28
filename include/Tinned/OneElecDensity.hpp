/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of one-electron spin-orbital density matrix.

   2023-09-21, Bin Gao:
   * first version
*/

#pragma once

#include "Tinned/ElectronicState.hpp"

namespace Tinned
{
    class OneElecDensity: public ElectronicState
    {
        public:
            //! Constructor
            // `derivatives` may be used only for `diff_impl()`
            explicit OneElecDensity(
                const std::string& name,
                const SymEngine::multiset_basic& derivatives = {}
            );

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;
    };

    // Helper function to make a one-electron spin-orbital density matrix
    inline SymEngine::RCP<const OneElecDensity> make_1el_density(
        const std::string& name
    )
    {
        return SymEngine::make_rcp<const OneElecDensity>(name);
    }
}
