/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of perturbation multichains.

   2024-10-31, Bin Gao:
   * changed to PertMultichain

   2024-04-30, Bin Gao:
   * first version
*/

#pragma once

#include <iostream>
#include "Tinned/StringifyVisitor.hpp"

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/symbol.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_rcp.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/ZeroOperator.hpp"
#include "Tinned/ZerosRemover.hpp"

namespace Tinned
{
    // Perturbation multichain
    typedef std::multiset<SymEngine::RCP<const Perturbation>, SymEngine::RCPBasicKeyLess>
        PertMultichain;

    // Convert derivatives to `PertMultichain`
    inline PertMultichain make_pert_multichain(
        const SymEngine::multiset_basic& derivatives
    )
    {
        PertMultichain perturbations;
        for (const auto& p: derivatives) {
            SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(*p))
            perturbations.insert(SymEngine::rcp_dynamic_cast<const Perturbation>(p));
        }
        return perturbations;
    }

    // Helper function to do high-order differentiation, and to remove zero
    // quantities
    template<typename T,
             typename std::enable_if<std::is_same<T, PertMultichain>::value ||
                 std::is_same<T, SymEngine::multiset_basic>::value ||
                 std::is_same<T, SymEngine::vec_basic>::value, int>::type = 0>
    inline SymEngine::RCP<const SymEngine::Basic> differentiate(
        const SymEngine::RCP<const SymEngine::Basic>& expr, const T& perturbations
    )
    {
        auto result = expr;
        for (const auto& p: perturbations) {
            SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(*p))
            result = result->diff(SymEngine::rcp_dynamic_cast<const Perturbation>(p));
        }
        result = remove_zeros(result);
        if (result.is_null()) {
            if (SymEngine::is_a_sub<const SymEngine::MatrixExpr>(*expr)) {
                return make_zero_operator();
            }
            else {
                return SymEngine::zero;
            }
        }
        else {
            return result;
        }
    }

    // Compute the sum of perturbation frequencies
    template<typename T,
             typename std::enable_if<std::is_same<T, PertMultichain>::value ||
                 std::is_same<T, SymEngine::multiset_basic>::value ||
                 std::is_same<T, SymEngine::vec_basic>::value, int>::type = 0>
    inline SymEngine::RCP<const SymEngine::Basic> get_frequency_sum(
        const T& perturbations
    )
    {
        SymEngine::RCP<const SymEngine::Basic> result = SymEngine::zero;
        for (const auto& p: perturbations) {
            SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(*p))
            result = SymEngine::add(
                result,
                SymEngine::rcp_dynamic_cast<const Perturbation>(p)->get_frequency()
            );
        }
        return result;
    }
}
