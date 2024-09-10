/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of time differentiation operator
   i\frac{\partial}{\partial t}.

   2024-06-04, Bin Gao:
   * remove the support for `NonElecFunction`

   2023-10-01, Bin Gao:
   * first version
*/

#pragma once

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/symengine_exception.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/matrix_symbol.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/OneElecOperator.hpp"

namespace Tinned
{
    // <i\frac{\partial}{\partial t}|, |i\frac{\partial}{\partial t}>
    enum TemporumType {Ket = 0, Bra = 1};

    class TemporumOperator: public SymEngine::MatrixSymbol
    {
        protected:
            // Target that the time differentiation operator acts on
            SymEngine::RCP<const SymEngine::MatrixExpr> target_;
            // Type of the time differentiation operator
            TemporumType type_;

        public:
            explicit TemporumOperator(
                const SymEngine::RCP<const SymEngine::MatrixExpr>& target,
                const TemporumType type
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get target
            inline SymEngine::RCP<const SymEngine::MatrixExpr> get_target() const
            {
                return target_;
            }

            // Get type of the time differentiation operator
            inline TemporumType get_type() const
            {
                return type_;
            }

            // Get frequency factor +/-sum(w)
            inline SymEngine::RCP<const SymEngine::Basic> get_frequency() const
            {
                SymEngine::RCP<const SymEngine::Basic> result = SymEngine::zero;
                for (const auto& var: get_derivatives()) {
                    auto pert = SymEngine::rcp_dynamic_cast<const Perturbation>(var);
                    result = SymEngine::add(result, pert->get_frequency());
                }
                return type_==TemporumType::Ket
                    ? result : SymEngine::sub(SymEngine::zero, result);
            }

            // Get derivatives
            inline SymEngine::multiset_basic get_derivatives() const
            {
                if (SymEngine::is_a_sub<const ElectronicState>(*target_)) {
                    auto target = SymEngine::rcp_dynamic_cast<const ElectronicState>(target_);
                    return target->get_derivatives();
                }
                else if (SymEngine::is_a_sub<const OneElecOperator>(*target_)) {
                    auto target = SymEngine::rcp_dynamic_cast<const OneElecOperator>(target_);
                    return target->get_derivatives();
                }
                else {
                    throw SymEngine::SymEngineException(
                        "Invalid type from the time differentiated target."
                    );
                }
            }
    };

    // Helper function to make a time differentiation operator
    inline SymEngine::RCP<const TemporumOperator> make_dt_operator(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& target,
        const TemporumType type = TemporumType::Ket
    )
    {
        return SymEngine::make_rcp<const TemporumOperator>(target, type);
    }
}
