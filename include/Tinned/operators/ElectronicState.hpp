/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of abstract class for electronic state.

   2023-10-28, Bin Gao:
   * remove member method get_args()

   2023-09-21, Bin Gao:
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
    // ElectronicState can be differentiated to any perturbation and any order
    class ElectronicState: public SymEngine::MatrixSymbol
    {
        protected:
            // derivatives_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivatives_;

        public:
            //! Constructor
            explicit ElectronicState(
                const std::string& name,
                const SymEngine::multiset_basic& derivatives = {}
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            //SymEngine::vec_basic get_args() const override;

            //// Override the defaut behaviour for diff
            //SymEngine::RCP<const SymEngine::Basic> diff_impl(
            //    const SymEngine::RCP<const SymEngine::Symbol>& s
            //) const override;

            // Get derivatives
            inline SymEngine::multiset_basic get_derivatives() const
            {
                return derivatives_;
            }

            // Check if `x` is a same response parameter
            virtual inline bool is_same_parameter(
                const SymEngine::RCP<const SymEngine::Basic>& x
            ) const = 0;
    };
}
