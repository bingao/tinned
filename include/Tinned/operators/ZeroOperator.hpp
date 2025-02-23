/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of zero operator, usually from differentiation.

   2024-06-05, Bin Gao:
   * inherits from `SymEngine::ZeroMatrix` so that `ZeroOperator` can be
     "recognized" by SymEngine

   2024-05-07, Bin Gao:
   * first version
*/

#pragma once

#include <string>

#include <symengine/basic.h>
#include <symengine/constants.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/zero_matrix.h>

namespace Tinned
{
    class ZeroOperator: public SymEngine::ZeroMatrix
    {
        protected:
            std::string name_;

        public:
            //! Constructor
            explicit ZeroOperator():
                // Dimensions are meaningless here
                SymEngine::ZeroMatrix(SymEngine::one, SymEngine::one),
                name_(std::string("0"))
            {
                SYMENGINE_ASSIGN_TYPEID()
            }

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;

            // Get name
            inline std::string get_name() const { return name_; }
    };

    // Helper function to make a zero operator
    inline SymEngine::RCP<const ZeroOperator> make_zero_operator()
    {
        return SymEngine::make_rcp<const ZeroOperator>();
    }
}
