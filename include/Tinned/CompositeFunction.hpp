/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of composite functions.

   2023-10-23, Bin Gao:
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

namespace Tinned
{
    class CompositeFunction: public SymEngine::FunctionWrapper
    {
        protected:
            // Order of differentiation of the outer function
            unsigned int order_;
            // Inner function
            SymEngine::RCP<const SymEngine::Basic> inner_;

        public:
            //! Constructor
            explicit CompositeFunction(
                const std::string& name,
                const SymEngine::RCP<const SymEngine::Basic> inner,
                const unsigned int order = 0
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
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            inline unsigned int get_order() const {
                return order_;
            }

            inline SymEngine::RCP<const SymEngine::Basic> get_inner() const {
                return inner_;
            }
    };
}
