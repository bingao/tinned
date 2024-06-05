/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This header file implements utilities that are particularly useful for
   response theory.

   2024-05-07, Bin Gao:
   * first version
*/

#pragma once

#include <map>

#include <symengine/basic.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/ZeroOperator.hpp"

#include "Tinned/FindAllVisitor.hpp"
#include "Tinned/ZerosRemover.hpp"
#include "Tinned/ReplaceVisitor.hpp"
#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    // Helper function to do high-order differentiation, and to remove zero
    // quantities
    inline SymEngine::RCP<const SymEngine::Basic> differentiate(
        const SymEngine::RCP<const SymEngine::Basic>& expr,
        const PerturbationTuple& perturbations
    )
    {
        auto result = expr;
        for (const auto& p: perturbations) result = result->diff(p);
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

    // Map for the substitution of Tinned objects (type `T`) with SymEngine
    // `Basic` symbols as well as their derivatives. The derivatives can be got
    // from the function `get_derivatives()` of those Tinned objects.
    template<typename T>
    using TinnedBasicMap = std::map<SymEngine::RCP<const T>,
                                    SymEngine::RCP<const SymEngine::Basic>,
                                    SymEngine::RCPBasicKeyLess>;

    // Helper function to replace Tinned objects and their derivatives with
    // SymEngine `Basic` symbols and corresponding derivatives
    template<typename T>
    inline SymEngine::RCP<const SymEngine::Basic> replace_all(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const TinnedBasicMap<T>& subs_dict
    )
    {
        SymEngine::map_basic_basic diff_subs_dict;
        // For each Tinned object `obj`, find all its derivatives in `x` and
        // that will be replaced with the corresponding symbol and its
        // derivatives
        for (const auto& d: subs_dict) {
            for (const auto& obj: find_all(x, d.first)) {
                auto diff_symbol = d.second;
                //FIXME: should `get_derivatives` return perturbations of type SymEngine::RCP<const Perturbation>?
                if (SymEngine::is_a_sub<const T>(*obj)) {
                    auto op = SymEngine::rcp_dynamic_cast<const T>(obj);
                    for (const auto& p: op->get_derivatives()) {
                        if (SymEngine::is_a_sub<const SymEngine::Symbol>(*p)) {
                            auto s = SymEngine::rcp_dynamic_cast<const SymEngine::Symbol>(p);
                            diff_symbol = diff_symbol->diff(s);
                        }
                        else {
                            throw SymEngine::SymEngineException(
                                "replace_all() gets an invalid perturbation "+stringify(p)
                            );
                        }
                    }
                }
                else {
                    throw SymEngine::SymEngineException(
                        "replace_all() gets an invalid object "+stringify(obj)
                    );
                }
                diff_subs_dict.insert({obj, diff_symbol});
            }
        }
        if (diff_subs_dict.empty()) return x;
        return replace(x, diff_subs_dict);
    }
}
