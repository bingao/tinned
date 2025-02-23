/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of one-electron like operators.

   2023-10-27, Bin Gao:
   * remove member method get_args()

   2023-09-08, Bin Gao:
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

namespace Tinned
{
    class OneElecOperator: public SymEngine::MatrixSymbol
    {
        protected:
            // dependencies_ stores perturbations that the operator depends on
            // and their maximum orders that can be differentiated
            PertDependency dependencies_;
            // derivatives_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivatives_;

        public:
            //! Constructor
            // `derivatives` may be used only for `diff_impl()`
            explicit OneElecOperator(
                const std::string& name,
                const PertDependency& dependencies,
                const SymEngine::multiset_basic& derivatives = {}
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            //SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

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

    // Helper function to make a one-electron like operator
    inline SymEngine::RCP<const OneElecOperator> make_1el_operator(
        const std::string& name,
        const PertDependency& dependencies = {}
    )
    {
        return SymEngine::make_rcp<const OneElecOperator>(name, dependencies);
    }
}
