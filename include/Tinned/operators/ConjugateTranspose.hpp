/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of conjugate tranpose.

   2024-06-04, Bin Gao:
   * first version
*/

#pragma once

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/identity_matrix.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/matrices/diagonal_matrix.h>
#include <symengine/matrices/immutable_dense_matrix.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/hadamard_product.h>
#include <symengine/visitor.h>

namespace Tinned
{
    // We have to modify SymEngine if ConjugateTranspose inherits from `MatrixExpr`
    class ConjugateTranspose: public SymEngine::MatrixSymbol
    {
        protected:
            SymEngine::RCP<const SymEngine::MatrixExpr> arg_;

        public:
            explicit ConjugateTranspose(
                const SymEngine::RCP<const SymEngine::MatrixExpr>& arg
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get the argument
            inline SymEngine::RCP<const SymEngine::MatrixExpr> get_arg() const
            {
                return arg_;
            }
    };

    // Mostly for internal use
    class ConjugateTransposeVisitor: public SymEngine::BaseVisitor<ConjugateTransposeVisitor>
    {
        private:
            SymEngine::RCP<const SymEngine::MatrixExpr> result_;

        public:
            explicit ConjugateTransposeVisitor() = default;

            void bvisit(const SymEngine::Basic& x);
            void bvisit(const SymEngine::MatrixExpr& x);
            void bvisit(const SymEngine::IdentityMatrix& x);
            void bvisit(const SymEngine::ZeroMatrix& x);
            void bvisit(const SymEngine::DiagonalMatrix& x);
            void bvisit(const SymEngine::ImmutableDenseMatrix& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixAdd& x);
            void bvisit(const SymEngine::HadamardProduct& x);

            inline SymEngine::RCP<const SymEngine::MatrixExpr> apply(
                const SymEngine::RCP<const SymEngine::MatrixExpr>& s
            )
            {
                s->accept(*this);
                return result_;
            }
    };

    // Helper function to make a conjugate transpose
    inline SymEngine::RCP<const SymEngine::MatrixExpr> make_conjugate_transpose(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& arg
    )
    {
        ConjugateTransposeVisitor visitor;
        return visitor.apply(arg);
    }
}
