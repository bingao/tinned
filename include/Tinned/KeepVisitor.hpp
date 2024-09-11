/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of keeping specific symbols while removing
   others.

   2024-05-10, Bin Gao:
   * add option to remove zero quantities in the function `keep_if`

   2023-10-19, Bin Gao:
   * moved from file Tinned/RemoveVisitor.hpp
*/

#pragma once

#include <functional>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/mul.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

#include "Tinned/ZerosRemover.hpp"
#include "Tinned/RemoveVisitor.hpp"

namespace Tinned
{
    // Keeping symbols if they match any given ones while removing others
    class KeepVisitor: public SymEngine::BaseVisitor<KeepVisitor, RemoveVisitor>
    {
        protected:
            // Check inequality for `x` and symbols to be kept
            inline bool is_not_equal(const SymEngine::Basic& x) const
            {
                for (const auto& s: symbols_) {
                    if (SymEngine::eq(x, *s)) return false;
                }
                return true;
            }

            // Function template for only one argument.
            template <typename Arg>
            inline void keep_if_arguments(
                SymEngine::vec_basic &f_args,
                bool &has_arg_kept,
                bool &has_arg_affected,
                const Arg &arg);

            // Function template for one or more arguments. `f_args` holds all
            // arguments, either affected or unaffected after removal.
            // `has_arg_kept` indicates if one or more arguments are kept.
            // `has_arg_affected` indicates if one or more arguments are
            // affected due to removal.
            template <typename FirstArg, typename... Args>
            inline void keep_if_arguments(
                SymEngine::vec_basic &f_args,
                bool &has_arg_kept,
                bool &has_arg_affected,
                const FirstArg &first_arg, const Args &...args)
            {
                keep_if_arguments(f_args, has_arg_kept, has_arg_affected, first_arg);
                keep_if_arguments(f_args, has_arg_kept, has_arg_affected, args...);
            }

            // Function template for a function like object with one or more arguments
            template<typename Fun, typename FirstArg, typename... Args>
            inline void keep_if_a_function(
                const Fun& x,
                const std::function<SymEngine::RCP<const SymEngine::Basic>(
                    const SymEngine::vec_basic&
                )>& constructor,
                const FirstArg& first_arg,
                const Args&... args
            )
            {
                // If the function will not be kept as whole, we then check if
                // any of its argument will be kept
                if (condition_(x)) {
                    auto f_args = SymEngine::vec_basic({});
                    auto has_arg_kept = false;
                    auto has_arg_affected = false;
                    keep_if_arguments(
                        f_args, has_arg_kept, has_arg_affected, first_arg, args...
                    );
                    if (has_arg_kept) {
                        result_ = has_arg_affected
                                ? constructor(f_args) : x.rcp_from_this();
                    }
                    // `result_` will be null only if all arguments are removed
                    else {
                        result_ = SymEngine::RCP<const SymEngine::Basic>();
                    }
                }
                // The function will be kept as a whole
                else {
                    result_ = x.rcp_from_this();
                }
            }

        public:
            explicit KeepVisitor(
                const SymEngine::set_basic& symbols
            ) : SymEngine::BaseVisitor<KeepVisitor, RemoveVisitor>(
                    symbols,
                    [&](const SymEngine::Basic& x) -> bool
                    {
                        return this->is_not_equal(x);
                    }
                )
            {
            }

            inline SymEngine::RCP<const SymEngine::Basic> apply(
                const SymEngine::RCP<const SymEngine::Basic>& x
            )
            {
                if (condition_(*x)) {
                    x->accept(*this);
                } else {
                    result_ = x;
                }
                return result_;
            }

            using RemoveVisitor::bvisit;
            //
            // Different from `RemoveVisitor`, the whole `Mul`, `MatrixMul`
            // and `HadamardProduct` will be kept whenever there is one factor
            // matches given symbols. Moreover, a function or an operator will
            // be kept if one of its argument matches given symbols.
            //
            void bvisit(const SymEngine::Add& x);
            void bvisit(const SymEngine::Mul& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            void bvisit(const SymEngine::Trace& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixAdd& x);
            void bvisit(const SymEngine::MatrixMul& x);
    };

    // Helper function to keep given `symbols` in `x` while removing others.
    // Note that zero quantities may produce after processing `MatrixMul`. One
    // can call `remove_zeros` on the result from `keep_if` if there are no
    // zero quantities in `x`.
    inline SymEngine::RCP<const SymEngine::Basic> keep_if(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::set_basic& symbols,
        const bool remove_zero_quantities = true
    )
    {
        KeepVisitor visitor(symbols);
        auto result = visitor.apply(x);
        if (result.is_null()) return result;
        return remove_zero_quantities ? remove_zeros(result) : result;
    }
    // Function template for only one argument.
    template <typename Arg>
    inline void KeepVisitor::keep_if_arguments(
        SymEngine::vec_basic &f_args,
        bool &has_arg_kept,
        bool &has_arg_affected,
        const Arg &arg)
    {
        auto new_arg = apply(arg);
        // If there exists an argument being kept, all other arguments
        // and the function will be kept. So we save all arguments to
        // be removed for later use.
        if (new_arg.is_null())
        {
            f_args.push_back(arg);
        }
        else
        {
            f_args.push_back(new_arg);
            has_arg_kept = true;
            if (SymEngine::neq(*arg, *new_arg))
                has_arg_affected = true;
        }
    }

    // Function template for only one argument of type `SymEngine::vec_basic`.
    template <>
    inline void KeepVisitor::keep_if_arguments<SymEngine::vec_basic>(
        SymEngine::vec_basic &f_args,
        bool &has_arg_kept,
        bool &has_arg_affected,
        const SymEngine::vec_basic &arg)
    {
        for (const auto &term : arg)
            keep_if_arguments(f_args, has_arg_kept, has_arg_affected, term);
    }
}
