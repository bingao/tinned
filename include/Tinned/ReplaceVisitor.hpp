/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of symbolic replacement.

   2024-04-29, Bin Gao:
   * add function `replace_with_derivatives()`

   2023-09-24, Bin Gao:
   * first version
*/

#pragma once

#include <map>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>
#include <symengine/subs.h>

#include "Tinned/FindAllVisitor.hpp"
#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    class ReplaceVisitor: public SymEngine::BaseVisitor<ReplaceVisitor, SymEngine::MSubsVisitor>
    {
        protected:
            // Template method that replaces `x` as a whole
            template<typename T> inline void replace_whole(T& x)
            {
                for (const auto& p: subs_dict_) {
                    if (SymEngine::eq(x, *p.first)) {
                        result_ = p.second;
                        return;
                    }
                }
                result_ = x.rcp_from_this();
            }

            // Template method for one argument function like classes
            template<typename Fun, typename Arg>
            inline void replace_one_arg_f(
                Fun& x,
                const SymEngine::RCP<Arg>& arg,
                std::function<SymEngine::RCP<const SymEngine::Basic>(
                    const SymEngine::RCP<Arg>&
                )> constructor
            )
            {
                // We first check if the argument will be replaced
                auto new_arg = apply(arg);
                SymEngine::RCP<const SymEngine::Basic> new_fun;
                if (SymEngine::eq(*arg, *new_arg)) {
                    new_fun = x.rcp_from_this();
                }
                else {
                    new_fun = constructor(
                        SymEngine::rcp_dynamic_cast<Arg>(new_arg)
                    );
                }
                // Next we check if the "new" function will be replaced
                replace_whole<Fun>(SymEngine::down_cast<const Fun&>(*new_fun));
            }

        public:
            explicit ReplaceVisitor(
                const SymEngine::map_basic_basic& subs_dict_,
                bool cache = true
            ) : SymEngine::BaseVisitor<ReplaceVisitor, SymEngine::MSubsVisitor>(
                    subs_dict_, cache
                )
            {
            }

            using SymEngine::MSubsVisitor::bvisit;
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::ZeroMatrix& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            void bvisit(const SymEngine::Trace& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixAdd& x);
            void bvisit(const SymEngine::MatrixMul& x);
            void bvisit(const SymEngine::MatrixDerivative& x);
    };

    // Helper function to replace classes defined in Tinned library in addition
    // to those in SymEngine::msubs()
    inline SymEngine::RCP<const SymEngine::Basic> replace(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::map_basic_basic& subs_dict
        //bool cache = false
    )
    {
        //ReplaceVisitor visitor(subs_dict, cache);
        ReplaceVisitor visitor(subs_dict, false);
        return visitor.apply(x);
    }

    // Map for the substitution of Tinned classes with SymEngine `Basic`
    // symbols as well as their derivatives
    template<typename T>
    using TinnedBasicMap = std::map<SymEngine::RCP<const T>,
                                    SymEngine::RCP<const SymEngine::Basic>,
                                    SymEngine::RCPBasicKeyLess>;

    // Helper function to replace classes defined in Tinned library as well as
    // their derivatives with symbols that they represent and corresponding
    // derivatives
    template<typename T>
    inline SymEngine::RCP<const SymEngine::Basic> replace_with_derivatives(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const TinnedBasicMap<T>& subs_dict
    )
    {
        SymEngine::map_basic_basic diff_subs_dict;
        // For each Tinned class `rep`, find all its derivatives in `x` and
        // that will be replaced with a symbol and its derivatives
        for (const auto& d: subs_dict) {
            for (const auto& rep: find_all<T>(x, d.first)) {
                auto diff_symbol = d.second;
                //FIXME: we may store derivatives in terms of SymEngine::RCP<const Perturbation>?
                for (const auto& p: rep->get_derivatives()) {
                    if (SymEngine::is_a_sub<const SymEngine::Symbol>(*p)) {
                        auto s = SymEngine::rcp_dynamic_cast<const SymEngine::Symbol>(p);
                        diff_symbol = diff_symbol->diff(s);
                    }
                    else {
                        throw SymEngine::SymEngineException(
                            "replace_with_derivatives() gets an invalid perturbation "
                            + stringify(p)
                        );
                    }
                }
                diff_subs_dict.insert({rep, diff_symbol});
            }
        }
        ReplaceVisitor visitor(diff_subs_dict, false);
        return visitor.apply(x);
    }
}
