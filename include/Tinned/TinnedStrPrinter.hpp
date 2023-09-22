/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of string printer.

   2023-09-19, Bin Gao:
   * first version
*/

#pragma once

#include <ostream>
#include <string>

#include <symengine/basic.h>
//#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>
#include <symengine/printers/strprinter.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecOperator.hpp"

namespace Tinned
{
    class TinnedStrPrinter: public SymEngine::BaseVisitor<TinnedStrPrinter, SymEngine::StrPrinter>
    {
        public:
            using SymEngine::StrPrinter::bvisit;
            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const FunctionSymbol& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            //std::string apply(const SymEngine::Basic& x);
    };

    // Stringify symbols from SymEngine and additional ones defined in this
    // library
    std::string tinned_str(const SymEngine::Basic& x);
    //std::string tinned_str(SymEngine::RCP<const SymEngine::Basic>& x);
}
