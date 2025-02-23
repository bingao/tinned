/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of adjoint maps in which all X's and all their
   derivatives are commutative.

   2024-04-17, Bin Gao:
   * first version
*/

#pragma once

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_symbol.h>

namespace Tinned
{
    class AdjointMap: public SymEngine::MatrixSymbol
    {
        protected:
            //FIXME: types of x_ and y_ should be `SymEngine::MatrixExpr`
            SymEngine::vec_basic x_;
            SymEngine::RCP<const SymEngine::Basic> y_;

        public:
            // One can use nested `AdjointMap`'s if X's are non-commutative
            explicit AdjointMap(
                const SymEngine::vec_basic& x,
                const SymEngine::RCP<const SymEngine::Basic>& y
            );

            // Constructor for `diff_impl()` of `ClusterConjHamiltonian`, where
            // `x` and its derivatives commute with X's and their derivatives
            // in `other`
            explicit AdjointMap(
                const AdjointMap& other,
                const SymEngine::RCP<const SymEngine::Basic>& x
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get the vector of operators X's
            inline SymEngine::vec_basic get_x() const
            {
                return x_;
            }

            // Get the size of X's
            inline std::size_t size_x() const
            {
                return x_.size();
            }

            // Get the operator Y
            inline SymEngine::RCP<const SymEngine::Basic> get_y() const
            {
                return y_;
            }
    };

    // Helper function to make an adjoint map
    inline SymEngine::RCP<const AdjointMap> make_adjoint_map(
        const SymEngine::vec_basic& x,
        const SymEngine::RCP<const SymEngine::Basic>& y
    )
    {
        return SymEngine::make_rcp<const AdjointMap>(x, y);
    }
}
