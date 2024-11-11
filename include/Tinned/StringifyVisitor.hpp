/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of string printer.

   2024-05-06, Bin Gao:
   * add verbose and compact modes

   2023-09-19, Bin Gao:
   * first version
*/

#pragma once

#include <ostream>
#include <regex>
#include <string>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/functions.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>
#include <symengine/printers/strprinter.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/ElectronicState.hpp"

#include "Tinned/FindAllVisitor.hpp"

namespace Tinned
{
    class StringifyVisitor: public SymEngine::BaseVisitor<StringifyVisitor, SymEngine::StrPrinter>
    {
        protected:
            // In verbose mode, functions/operators will be printed with their
            // dependencies, perturbations will be printed with their
            // components and frequencies
            bool verbose_;

            // Put `expr` into square brackets
            inline std::string square_bracket(const std::string& expr) const
            {
                return "[" + expr + "]";
            }

            // Put `expr` into angle brackets
            inline std::string angle_bracket(const std::string& expr) const
            {
                return "<" + expr + ">";
            }

            // Stringify an operator
            inline std::string stringify_operator(
                const std::string& name,
                const SymEngine::multiset_basic& derivatives = {},
                const PertDependency& dependencies = {}
            )
            {
                std::ostringstream o;
                if (derivatives.empty()) {
                    if (dependencies.empty() || !verbose_) return name;
                    for (auto p=dependencies.begin(); p!=dependencies.end(); ++p) {
                        if (p==dependencies.begin()) {
                            o << apply(*p->first)
                              << angle_bracket(std::to_string(p->second));
                        }
                        else {
                            o << ", " << apply(*p->first)
                              << angle_bracket(std::to_string(p->second));
                        }
                    }
                    return name + parenthesize(o.str());
                }
                else {
                    o << name;
                    for (const auto& var: derivatives) o << ", " << apply(var);
                    return "Derivative" + parenthesize(o.str());
                }
            }

            // Stringify an electronic state
            inline std::string stringify_state(
                const SymEngine::RCP<const ElectronicState>& state
            )
            {
                bvisit(*state);
                return str_;
            }

        public:
            explicit StringifyVisitor(const bool verbose = true): verbose_(verbose) {}

            using SymEngine::StrPrinter::bvisit;
            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::ZeroMatrix& x);
            void bvisit(const SymEngine::MatrixSymbol& x);

            // Convert `x` into a string
            inline std::string convert(const SymEngine::Basic& x)
            {
                x.accept(*this);
                // Replace "+-1*" with "-", and "+-" with "-"
                str_ = std::regex_replace(str_, std::regex(R"(\+-1\*)"), R"(-)");
                str_ = std::regex_replace(str_, std::regex(R"(\+-)"), R"(-)");
                return str_;
            }
    };

    // Helper functions to stringify symbols from SymEngine and additional ones
    // defined in Tinned library
    inline std::string stringify(const SymEngine::Basic& x, const bool verbose = true)
    {
        StringifyVisitor visitor(verbose);
        return visitor.convert(x);
    }

    inline std::string stringify(
        const SymEngine::RCP<const SymEngine::Basic>& x, const bool verbose = true
    )
    {
        StringifyVisitor visitor(verbose);
        return visitor.convert(*x);
    }

    // Stringify the result from `FindAllVisitor`
    //FIXME: incorrect new line character
    inline std::string stringify(const FindAllResult& result)
    {
        auto str_result = std::string("{\\\n");
        for (const auto& term: result) {
            str_result += "    {" + std::to_string(term.first) + ", {";
            for (const auto& op: term.second) str_result += ", " + stringify(op);
            str_result += "}\\\n    },\\\n";
        }
        return str_result + "}";
    }
}
