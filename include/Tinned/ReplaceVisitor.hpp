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
#include <symengine/subs.h>

namespace Tinned
{
    class ReplaceVisitor: public SymEngine::BaseVisitor<ReplaceVisitor, SymEngine::MSubsVisitor>
    {
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
            void bvisit(const SymEngine::Symbol& x);
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

    // Replace classes defined in this library in addition to those in
    // SymEngine::msubs()
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
