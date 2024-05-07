/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of zero operator, usually from differentiation.

   2024-05-07, Bin Gao:
   * first version
*/

#pragma once

#include <string>

#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_symbol.h>

//#include "Tinned/PertDependency.hpp"

namespace Tinned
{
    class ZeroOperator: public SymEngine::MatrixSymbol
    {
        public:
            //! Constructor
            explicit ZeroOperator(): SymEngine::MatrixSymbol(std::string("0"))
            {
                SYMENGINE_ASSIGN_TYPEID()
            }

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            //// Get dependencies
            //inline PertDependency get_dependencies() const { return {}; }

            //// Get derivatives
            //inline SymEngine::multiset_basic get_derivatives() const { return {}; }
    };

    // Helper function to make a zero operator
    inline SymEngine::RCP<const ZeroOperator> make_zero_operator()
    {
        return SymEngine::make_rcp<const ZeroOperator>();
    }
}
