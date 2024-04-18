/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

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
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>
#include <symengine/subs.h>

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
                std::function<SymEngine::RCP<Fun>(const SymEngine::RCP<Arg>&)> constructor
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
}
