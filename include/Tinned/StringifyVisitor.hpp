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
#include <symengine/dict.h>
//#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>
#include <symengine/printers/strprinter.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

namespace Tinned
{
    class StringifyVisitor: public SymEngine::BaseVisitor<StringifyVisitor, SymEngine::StrPrinter>
    {
        protected:
            // Stringify perturbation dependencies
            inline std::string to_string(
                const std::string& name,
                const PertDependency& dependencies
            )
            {
                std::ostringstream o;
                o << name << "(";
                for (auto p = dependencies.begin(); p != dependencies.end(); ++p) {
                    if (p == dependencies.begin()) {
                        o << apply(*p->first) << "<" << p->second << ">";
                    }
                    else {
                        o << ", " << apply(*p->first) << "<" << p->second << ">";
                    }
                }
                o << ")";
                return o.str();
            }

            // Stringify the derivative of an operator
            inline std::string to_string(
                const std::string& name,
                const SymEngine::multiset_basic& derivative
            )
            {
                std::ostringstream o;
                o << "Derivative(" << name;
                for (const auto& var: derivative) {
                    o << ", " << apply(*var);
                }
                o << ")";
                return o.str();
            }

            // Stringify the electronic state
            inline std::string to_string(
                const SymEngine::RCP<const ElectronicState>& state
            )
            {
                bvisit(*state);
                return str_;
            }

            // Stringify the derivative of an operator with an electronic state
            // as the argument
            inline std::string to_string(
                const std::string& name,
                const SymEngine::RCP<const ElectronicState>& state,
                const SymEngine::multiset_basic& derivative
            )
            {
                auto str_state = to_string(state);
                auto str_op = name + "(" + str_state + ")";
                return derivative.empty() ? str_op : to_string(str_op, derivative);
            }

        public:
            using SymEngine::StrPrinter::bvisit;
            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            //std::string apply(const SymEngine::Basic& x);
    };

    // Stringify symbols from SymEngine and additional ones defined in this
    // library
    std::string stringify(const SymEngine::Basic& x);
    //std::string stringify(SymEngine::RCP<const SymEngine::Basic>& x);
}
