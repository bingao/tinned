/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of (perturbed) response parameters that depend
   on all different perturbations at arbitrary order.

   2024-06-04, Bin Gao:
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
    class PerturbedParameter: public SymEngine::MatrixSymbol
    {
        protected:
            // derivatives_ holds derivatives with respect to perturbations
            SymEngine::multiset_basic derivatives_;

        public:
            //! Constructor
            // `derivatives` may be used only for `diff_impl()`
            explicit PerturbedParameter(
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
                if (SymEngine::is_a_sub<const PerturbedParameter>(*x)) {
                    auto op = SymEngine::rcp_dynamic_cast<const PerturbedParameter>(x);
                    return get_name()==op->get_name() ? true : false;
                }
                return false;
            }
    };

    // Helper function to make a (perturbed) response parameter
    inline SymEngine::RCP<const PerturbedParameter> make_perturbed_parameter(
        const std::string& name,
        const SymEngine::multiset_basic& derivatives = {}
    )
    {
        return SymEngine::make_rcp<const PerturbedParameter>(name, derivatives);
    }
}
