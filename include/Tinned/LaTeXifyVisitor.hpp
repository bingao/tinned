/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of LaTeX printer.

   2024-05-04, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>
#include <ostream>
#include <regex>
#include <string>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/functions.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>
#include <symengine/printers/latex.h>

#include "Tinned/ElectronicState.hpp"

namespace Tinned
{
    // Font style for operators
    enum OperFontStyle {Regular = 0, Bold = 1};

    class LaTeXifyVisitor: public SymEngine::BaseVisitor<LaTeXifyVisitor, SymEngine::LatexPrinter>
    {
        protected:
            // String for a newline
            const std::string str_newline_;
            // Maximum number of functions/operators per line
            unsigned int max_num_symbols_;
            // Current number of functions/operators per line
            unsigned int curr_num_symbols_;
            // Indicates if zero factor(s) exists in a multiplication
            bool zero_factor_;

            // Remove newline from `expr`
            inline std::string remove_newline(const std::string& expr) const
            {
                return std::regex_replace(expr, std::regex(str_newline_.data()), R"()");
            }

            // Find the position of trailing newline in `expr`
            inline std::size_t find_trailing_newline(const std::string& expr) const
            {
                auto pos = expr.rfind(str_newline_);
                return pos==expr.size()-str_newline_.size() ? pos : std::string::npos;
            }

            // Add a suffix to `expr`
            inline std::string add_suffix(
                const std::string& expr, const std::string& suffix
            ) const
            {
                auto pos = find_trailing_newline(expr);
                if (pos==std::string::npos) {
                    return expr + suffix;
                }
                else {
                    return expr.substr(0, pos) + suffix + str_newline_;
                }
            }

            // Update the current number of functions/operators per line and
            // make a newline if needed
            inline void update_num_symbols(
                const unsigned int num_symbols,
                std::string& expr
            )
            {
                curr_num_symbols_ += num_symbols;
                if (curr_num_symbols_>=max_num_symbols_) {
                    expr += str_newline_;
                    curr_num_symbols_ = 0;
                }
            }

            // LaTeXify derivatives
            inline std::string latexify_derivatives(
                const SymEngine::multiset_basic& derivatives
            )
            {
                // Important to use a new visitor to convert `derivatives`
                // because we don't want to change the number of symbols when
                // invoking `apply` method
                LaTeXifyVisitor visitor;
                std::ostringstream o;
                auto var = derivatives.begin();
                o << "^{" << visitor.apply(*var);
                ++var;
                for (; var!=derivatives.end(); ++var) o << "," << visitor.apply(*var);
                o << "}";
                return remove_newline(o.str());
            }

            // LaTeXify an operator with/without derivatives
            inline std::string latexify_operator(
                const std::string& name,
                const SymEngine::multiset_basic& derivatives = {},
                const OperFontStyle style = OperFontStyle::Bold
            )
            {
                std::ostringstream o;
                if (style==OperFontStyle::Bold) {
                    o << "\\bm{" << name << "}";
                }
                else {
                    o << name;
                }
                if (!derivatives.empty()) o << latexify_derivatives(derivatives);
                return o.str();
            }

            std::string parenthesize(const std::string& expr) override;
            // Add \times at the beginning of a newline for multiplications
            std::string print_mul() override;

        public:
            explicit LaTeXifyVisitor(const unsigned int max_num_symbols = 12)
              : str_newline_("\\\\\n&"),
                max_num_symbols_(max_num_symbols),
                curr_num_symbols_(0),
                zero_factor_(false) {}

            using SymEngine::LatexPrinter::bvisit;
            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::Add& x);
            void bvisit(const SymEngine::Mul& x);
            void bvisit(const SymEngine::ZeroMatrix& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            void bvisit(const SymEngine::MatrixAdd& x);
            void bvisit(const SymEngine::MatrixMul& x);
            void bvisit(const SymEngine::Trace& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixDerivative& x);

            // Convert `x` into a LaTeX string
            inline std::string convert(const SymEngine::Basic& x)
            {
                x.accept(*this);
                // Remove trailing newline and add `&` for multiple lines
                auto pos = find_trailing_newline(str_);
                if (pos!=std::string::npos) str_.erase(pos, str_newline_.size());
                return str_.find(str_newline_)==std::string::npos ? str_ : "&" + str_;
            }
    };

    // Helper functions to latexify symbols from SymEngine and additional ones
    // defined in Tinned library
    inline std::string latexify(
        const SymEngine::Basic& x,
        const unsigned int max_num_symbols = 12
    )
    {
        LaTeXifyVisitor visitor(max_num_symbols);
        return visitor.convert(x);
    }

    inline std::string latexify(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const unsigned int max_num_symbols = 12
    )
    {
        LaTeXifyVisitor visitor(max_num_symbols);
        return visitor.convert(*x);
    }
}
