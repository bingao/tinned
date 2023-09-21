/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of perturbation.

   2023-09-16, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>

#include <symengine/dict.h>
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

namespace Tinned
{
    class Perturbation: public SymEngine::Symbol
    {
        private:
            std::size_t dimension_;

        public:
            //! Constructor
            explicit Perturbation(const char name, const std::size_t dimension);

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;

            //SymEngine::vec_basic get_args() const override
            //{
            //    return SymEngine::vec_basic({dimension_});
            //}

            //! Get dimension of the perturbation
            inline std::size_t get_dimension() const
            {
                return dimension_;
            }
    };
}
