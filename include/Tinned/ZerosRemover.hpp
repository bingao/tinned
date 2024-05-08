/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of removing zero quantities.

   2024-05-08, Bin Gao:
   * first version
*/

#pragma once

#include <functional>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/dict.h>
#include <symengine/mul.h>
#include <symengine/number.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

#include "Tinned/ZeroOperator.hpp"

namespace Tinned
{
    // Helper function to check if a quantity is the type `ZeroOperator` or a
    // number equal to zero
    inline bool is_zero_quantity(const SymEngine::Basic& x)
    {
        if (SymEngine::is_a_sub<const ZeroOperator>(x) ||
            SymEngine::is_a_sub<const SymEngine::ZeroMatrix>(x)) {
            return true;
        }
        else {
            return SymEngine::is_number_and_zero(x);
        }
    }

    // Different from `RemoveVisitor`, we check if a quantity is zero by
    // calling `is_zero_quantity()` instead of asking users to provide
    // `ZeroOperator`, zero constants (integer, real, complex, etc.) and
    // matrices.
    class ZerosRemover: public SymEngine::BaseVisitor<ZerosRemover>
    {
        protected:
            SymEngine::RCP<const SymEngine::Basic> result_;

            // Template method for one argument function like classes
            template<typename Fun, typename Arg>
            inline void remove_one_arg_f(
                Fun& x,
                const SymEngine::RCP<Arg>& arg,
                std::function<SymEngine::RCP<const SymEngine::Basic>(
                    const SymEngine::RCP<Arg>&
                )> constructor
            )
            {
                auto new_arg = apply(arg);
                if (new_arg.is_null()) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    if (SymEngine::eq(*arg, *new_arg)) {
                        result_ = x.rcp_from_this();
                    }
                    else {
                        result_ = constructor(
                            SymEngine::rcp_dynamic_cast<Arg>(new_arg)
                        );
                    }
                }
            }

        public:
            explicit ZerosRemover() noexcept = default;

            inline SymEngine::RCP<const SymEngine::Basic> apply(
                const SymEngine::RCP<const SymEngine::Basic>& x
            )
            {
                x->accept(*this);
                return result_;
            }

            void bvisit(const SymEngine::Basic& x);
            void bvisit(const SymEngine::Add& x);
            void bvisit(const SymEngine::Mul& x);
            void bvisit(const SymEngine::Trace& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixAdd& x);
            void bvisit(const SymEngine::MatrixMul& x);
            void bvisit(const SymEngine::MatrixDerivative& x);
    };

    // Helper function to remove zero quantities from `x`
    inline SymEngine::RCP<const SymEngine::Basic> remove_zeros(
        const SymEngine::RCP<const SymEngine::Basic>& x
    )
    {
        ZerosRemover visitor;
        return visitor.apply(x);
    }
}
