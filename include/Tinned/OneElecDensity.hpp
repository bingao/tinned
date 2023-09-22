/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of one-electron spin-orbital density matrix.

   2023-09-21, Bin Gao:
   * first version
*/

#pragma once

#include "Tinned/ElectronState.hpp"

namespace Tinned
{
    class OneElecDensity: public ElectronState
    {
        public:
            //! Constructor
            explicit OneElecDensity(
                const std::string& name,
                const SymEngine::multiset_basic& derivative = {}
            );

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::MatrixExpr> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;
    };
}
