/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of symbolic replacement.

   2023-09-24, Bin Gao:
   * first version
*/

#pragma once

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>
#include <symengine/subs.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"

namespace Tinned
{
    class ReplaceVisitor: public SymEngine::BaseVisitor<ReplaceVisitor, SymEngine::XReplaceVisitor>
    {
        public:
            using SymEngine::XReplaceVisitor::bvisit;
            explicit ReplaceVisitor(
                const SymEngine::map_basic_basic& subs_dict_,
                bool cache = true
            ) : SymEngine::BaseVisitor<ReplaceVisitor, SymEngine::XReplaceVisitor>(
                    subs_dict_, cache
                )
            {
            }

            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
    };

    // Replace classes defined in this library in addition to those in
    // SymEngine::xreplace()
    inline SymEngine::RCP<const SymEngine::Basic> replace(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::map_basic_basic& subs_dict,
        bool cache
    )
    {
        ReplaceVisitor visitor(subs_dict, cache);
        return visitor.apply(x);
    }
}
