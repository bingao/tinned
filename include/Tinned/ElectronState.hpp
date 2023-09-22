/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of abstract class for electron state.

   2023-09-21, Bin Gao:
   * first version
*/

#pragma once

#include <string>

#include <symengine/dict.h>
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/matrix_symbol.h>

namespace Tinned
{
    // ElectronState can be differentiated to any perturbation and any order
    class ElectronState: public SymEngine::MatrixSymbol
    {
        protected:
            // derivative_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivative_;

        public:
            //! Constructor
            explicit ElectronState(
                const std::string& name,
                const SymEngine::multiset_basic& derivative = {}
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            //// Override the defaut behaviour for diff
            //SymEngine::RCP<const SymEngine::MatrixExpr> diff_impl(
            //    const SymEngine::RCP<const SymEngine::Symbol>& s
            //) const override;

            // Get derivative
            inline SymEngine::multiset_basic get_derivative() const
            {
                return derivative_;
            }
    };
}
