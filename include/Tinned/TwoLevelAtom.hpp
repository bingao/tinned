/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of two-level atom system.

   2024-07-15, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/number.h>
#include <symengine/integer.h>
#include <symengine/constants.h>
#include <symengine/derivative.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/matrices/immutable_dense_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"
#include "Tinned/OperatorEvaluator.hpp"
#include "Tinned/FunctionEvaluator.hpp"

// This two-level atom system is adopted from section 5.10.3, "Principles and
// Practices of Molecular Properties: Theory, Modeling and Simulations",
// Patrick Norman, Kenneth Ruud and Trond Saue.
namespace Tinned
{
    // Evaluator for different (electron) operators
    class TwoLevelOperator: public OperatorEvaluator<SymEngine::RCP<const SymEngine::MatrixExpr>>
    {
        protected:
            // Maximum allowed order
            unsigned int max_order_;
            // Unperturbed Hamiltonian and its value
            std::pair<SymEngine::RCP<const OneElecOperator>,
                      SymEngine::RCP<const SymEngine::MatrixExpr>> H0_;
            // Density matrix and the value of unperturbed one
            std::pair<SymEngine::RCP<const OneElecOperator>,
                      SymEngine::RCP<const SymEngine::MatrixExpr>> rho0_;
            // External field's operators and their values, each operator
            // should depend only on one unique perturbation
            std::map<SymEngine::RCP<const OneElecOperator>,
                     SymEngine::RCP<const SymEngine::MatrixExpr>,
                     SymEngine::RCPBasicKeyLess> V_;
            // Transition angular frequency matrix
            SymEngine::RCP<const SymEngine::ImmutableDenseMatrix> omega_;
            // All perturbations
            SymEngine::set_basic perturbations_;
            // Cached derivatives of density matrix <order, <perturbations, derivatives>>
            //FIXME: not optimal for identical perturbations
            typedef std::vector<std::pair<SymEngine::multiset_basic,
                                          SymEngine::RCP<const SymEngine::MatrixExpr>>>
                DensityDerivative;
            std::map<unsigned int, DensityDerivative> rho_all_derivatives_;

            // Make a 2x2 zero matrix
            inline SymEngine::RCP<const SymEngine::MatrixExpr> make_zero_matrix()
            {
                return SymEngine::immutable_dense_matrix(
                    2, 2, {SymEngine::zero, SymEngine::zero, SymEngine::zero, SymEngine::zero}
                );
            }

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_hermitian_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override;

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_1el_density(const OneElecDensity& x) override;

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_1el_operator(const OneElecOperator& x) override;

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_temporum_operator(const TemporumOperator& x) override;

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_conjugate_matrix(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override;

            SymEngine::RCP<const SymEngine::MatrixExpr>
            eval_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override;

            void eval_oper_addition(
                SymEngine::RCP<const SymEngine::MatrixExpr>& A,
                const SymEngine::RCP<const SymEngine::MatrixExpr>& B
            ) override;

            SymEngine::RCP<const SymEngine::MatrixExpr> eval_oper_multiplication(
                const SymEngine::RCP<const SymEngine::MatrixExpr>& A,
                const SymEngine::RCP<const SymEngine::MatrixExpr>& B
            ) override;

            void eval_oper_scale(
                const SymEngine::RCP<const SymEngine::Number>& scalar,
                SymEngine::RCP<const SymEngine::MatrixExpr>& A
            ) override;

        public:
            explicit TwoLevelOperator(
                const std::pair<SymEngine::RCP<const OneElecOperator>,
                                SymEngine::RCP<const SymEngine::MatrixExpr>>& H0,
                const std::map<SymEngine::RCP<const OneElecOperator>,
                               SymEngine::RCP<const SymEngine::MatrixExpr>,
                               SymEngine::RCPBasicKeyLess>& V,
                const std::pair<SymEngine::RCP<const OneElecOperator>,
                                SymEngine::RCP<const SymEngine::MatrixExpr>>& rho0
            );

            ~TwoLevelOperator() = default;
    };

    // Evaluator for different expectation values
    class TwoLevelFunction: public FunctionEvaluator<SymEngine::RCP<const SymEngine::Basic>,
                                                     SymEngine::RCP<const SymEngine::MatrixExpr>>
    {
        protected:
            SymEngine::RCP<const SymEngine::Basic> eval_trace(
                const SymEngine::RCP<const SymEngine::MatrixExpr>& A
            ) override;

            void eval_fun_addition(
                SymEngine::RCP<const SymEngine::Basic>& f,
                const SymEngine::RCP<const SymEngine::Basic>& g
            ) override;

            void eval_fun_scale(
                const SymEngine::RCP<const SymEngine::Number>& scalar,
                SymEngine::RCP<const SymEngine::Basic>& f
            ) override;

        public:
            explicit TwoLevelFunction(
                const std::pair<SymEngine::RCP<const OneElecOperator>,
                                SymEngine::RCP<const SymEngine::MatrixExpr>>& H0,
                const std::map<SymEngine::RCP<const OneElecOperator>,
                               SymEngine::RCP<const SymEngine::MatrixExpr>,
                               SymEngine::RCPBasicKeyLess>& V,
                const std::pair<SymEngine::RCP<const OneElecOperator>,
                                SymEngine::RCP<const SymEngine::MatrixExpr>>& rho0
            );

            ~TwoLevelFunction() = default;
    };
}
