/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of symbolic replacement.

   2024-06-13, Bin Gao:
   * move helper function template `replace_all` here

   2024-05-12, Bin Gao:
   * implement function templates `replace_a_function` and `replace_arguments`
     for the replacement of a function like object with one or more arguments

   2023-09-24, Bin Gao:
   * first version
*/

#pragma once

#include <map>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>
#include <symengine/subs.h>

#include "Tinned/PertMultichain.hpp"
#include "Tinned/FindAllVisitor.hpp"

namespace Tinned
{
    class ReplaceVisitor: public SymEngine::BaseVisitor<ReplaceVisitor, SymEngine::MSubsVisitor>
    {
        protected:
            // Function template that replaces `x` as a whole
            template<typename T> inline bool replace_a_whole(T& x)
            {
                for (const auto& p: subs_dict_) {
                    if (SymEngine::eq(x, *p.first)) {
                        result_ = p.second;
                        return true;
                    }
                }
                result_ = x.rcp_from_this();
                return false;
            }

            // Function template for replacing one argument.
            template <typename Arg>
            inline void replace_arguments(
                SymEngine::vec_basic &f_args,
                bool &has_arg_replaced,
                const Arg &arg)
            {
                auto new_arg = apply(arg);
                f_args.push_back(new_arg);
                if (SymEngine::neq(*arg, *new_arg))
                    has_arg_replaced = true;
            }

            // Function template for replacing one or more arguments. `f_args`
            // holds all arguments, either replaced or original ones.
            // `has_arg_replaced` indicates if one or more arguments are
            // replaced.
            template<typename FirstArg, typename... Args>
            inline void replace_arguments(
                SymEngine::vec_basic& f_args,
                bool& has_arg_replaced,
                const FirstArg& first_arg, const Args&... args
            )
            {
                replace_arguments(f_args, has_arg_replaced, first_arg);
                replace_arguments(f_args, has_arg_replaced, args...);
            }

            // Function template for a function like object with one or more arguments
            template<typename Fun, typename FirstArg, typename... Args>
            inline void replace_a_function(
                const Fun& x,
                const std::function<SymEngine::RCP<const SymEngine::Basic>(
                    const SymEngine::vec_basic&
                )>& constructor,
                const FirstArg& first_arg,
                const Args&... args
            )
            {
                // We first check if the function will be replaced as a whole
                if (!replace_a_whole(x)) {
                    // If the function is not replaced, we then perform the
                    // replacement of its arguments
                    auto f_args = SymEngine::vec_basic({});
                    auto has_arg_replaced = false;
                    replace_arguments(f_args, has_arg_replaced, first_arg, args...);
                    if (has_arg_replaced) {
                        result_ = constructor(f_args);
                    }
                    else {
                        result_ = x.rcp_from_this();
                    }
                }
            }

        public:
            explicit ReplaceVisitor(
                const SymEngine::map_basic_basic& subs_dict_,
                bool cache = false
            ) : SymEngine::BaseVisitor<ReplaceVisitor, SymEngine::MSubsVisitor>(
                    subs_dict_, cache
                )
            {}

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

    // Map for the substitution of Tinned objects (type `T`) with SymEngine
    // `Basic` symbols as well as their derivatives. The derivatives can be got
    // from the function `get_derivatives()` of those Tinned objects.
    template<typename T>
    using ReplaceDictMap = std::map<SymEngine::RCP<const T>,
                                    SymEngine::RCP<const SymEngine::Basic>,
                                    SymEngine::RCPBasicKeyLess>;

    // Helper function to replace Tinned objects and their derivatives with
    // SymEngine `Basic` symbols and corresponding derivatives
    template<typename T>
    inline SymEngine::RCP<const SymEngine::Basic> replace_all(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const ReplaceDictMap<T>& subs_dict
    )
    {
        SymEngine::map_basic_basic diff_subs_dict;
        // For each Tinned object `obj`, find all its derivatives in `x` and
        // that will be replaced with the corresponding symbol and its
        // derivatives
        for (const auto& d: subs_dict) {
            for (const auto& obj_map: find_all(x, d.first)) {
                for (const auto& obj: obj_map.second) {
                    SYMENGINE_ASSERT(SymEngine::is_a_sub<const T>(*obj))
                    auto op = SymEngine::rcp_dynamic_cast<const T>(obj);
                    diff_subs_dict.insert({
                        obj, differentiate(d.second, op->get_derivatives())
                    });
                }
            }
        }
        if (diff_subs_dict.empty()) return x;
        return replace(x, diff_subs_dict);
    }

    // Function template for only one argument of type `SymEngine::vec_basic`.
    template <>
    inline void ReplaceVisitor::replace_arguments<SymEngine::vec_basic>(
        SymEngine::vec_basic &f_args,
        bool &has_arg_replaced,
        const SymEngine::vec_basic &arg)
    {
        for (const auto &term : arg) replace_arguments(f_args, has_arg_replaced, term);
    }
}
