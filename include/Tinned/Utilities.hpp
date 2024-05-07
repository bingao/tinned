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
#include <symengine/matrices/matrix_mul.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/ZeroOperator.hpp"

#include "Tinned/FindAllVisitor.hpp"
#include "Tinned/RemoveVisitor.hpp"
#include "Tinned/ReplaceVisitor.hpp"
#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    // Helper function to do high-order differentiation, and to remove
    // `ZeroOperator` objects
    inline SymEngine::RCP<const SymEngine::Basic> differentiate(
        const SymEngine::RCP<const SymEngine::Basic>& expr,
        const PerturbationTuple& perturbations
    )
    {
        auto result = expr;
        for (const auto& p: perturbations) result = result->diff(p);
        return remove_if(
            result, SymEngine::set_basic({make_zero_operator(), SymEngine::zero})
        );
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
        ReplaceVisitor visitor(diff_subs_dict, false);
        return visitor.apply(x);
    }

    // Helper function to remove undifferentiated `TemporumOperator` objects
    // and replace differentiated ones with its corresponding target's
    // derivatives multiplied by sums of perturbation frequencies
    template<typename T>
    inline SymEngine::RCP<const SymEngine::Basic> replace_temporum(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::RCP<const T>& target,
        const TemporumType type = TemporumType::Ket
    )
    {
        auto dt_target = make_dt_operator(target, type);
        // Remove undifferentiated `TemporumOperator` objects in `x`
        auto diff_x = remove_if(x, SymEngine::set_basic({dt_target}));
        // Find all differentiated `TemporumOperator` objects in `diff_x`
        SymEngine::map_basic_basic dt_subs_dict;
        // `TemporumOperator` objects with zero sum of perturbation frequencies
        // should be removed before replacement
        SymEngine::set_basic dt_zero_freq;
        for (const auto& diff_dt: find_all(diff_x, dt_target)) {
            //FIXME: should `get_derivatives` return perturbations of type SymEngine::RCP<const Perturbation>?
            if (SymEngine::is_a_sub<const TemporumOperator>(*diff_dt)) {
                auto op = SymEngine::rcp_dynamic_cast<const TemporumOperator>(diff_dt);
                if (op->get_frequency()->is_zero()) {
                    dt_zero_freq.insert(diff_dt);
                }
                else {
                    SymEngine::RCP<const SymEngine::Basic> diff_target = target;
                    for (const auto& p: op->get_derivatives()) {
                        if (SymEngine::is_a_sub<const SymEngine::Symbol>(*p)) {
                            auto s = SymEngine::rcp_dynamic_cast<const SymEngine::Symbol>(p);
                            diff_target = diff_target->diff(s);
                        }
                        else {
                            throw SymEngine::SymEngineException(
                                "replace_temporum() gets an invalid perturbation "+stringify(p)
                            );
                        }
                    }
                    dt_subs_dict.insert({
                        diff_dt, SymEngine::matrix_mul({op->get_frequency(), diff_target})
                    });
                }
            }
            else {
                throw SymEngine::SymEngineException(
                    "replace_temporum() gets an invalid temporum "+stringify(diff_dt)
                );
            }
        }
        if (!dt_zero_freq.empty()) diff_x = remove_if(diff_x, dt_zero_freq);
        if (dt_subs_dict.empty()) return diff_x;
        ReplaceVisitor visitor(dt_subs_dict, false);
        return visitor.apply(diff_x);
//FIXME: after the above replacement, there might be new zero frequency terms due to the collection of same terms
    }
}
