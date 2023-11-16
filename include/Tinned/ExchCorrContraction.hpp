/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file implements functions for contractions between exchange-correlation
   energy functional derivative vectors and generalized density vectors.

   2023-11-08, Bin Gao:
   * first version
*/

#pragma once

#include <map>
#include <tuple>
#include <utility>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/symengine_rcp.h>
#include <symengine/symengine_assert.h>

#include "Tinned/NonElecFunction.hpp"
#include "Tinned/OneElecOperator.hpp"

namespace Tinned
{
    // Type for contractions between XC functional derivative vectors and
    // perturbed generalized density vectors
    typedef std::map<unsigned int, SymEngine::RCP<const SymEngine::Basic>>
        ExcDensityContractionMap;

    // Type for the collection of (un)perturbed weights, and corresponding
    // contractions between XC functional derivative vectors and perturbed
    // generalized density vectors
    //
    // Probably for developers: `term` is used as the iterator for the outer
    // map (`term.first` points to the grid weight and `term.second` to the
    // inner map), while `contr`  for the inner map (`contr.first` points to
    // the order of XC functional derivative and `contr.second` to the
    // perturbed generalized density vectors)
    typedef std::map<SymEngine::RCP<const NonElecFunction>, ExcDensityContractionMap,
                     SymEngine::RCPBasicKeyLess>
        ExcContractionMap;

    // Type for the collection of (un)perturbed weights, XC functional
    // derivative vectors, perturbed generalized density vectors and
    // (un)perturbed generalized overlap distributions
    typedef std::map<SymEngine::RCP<const OneElecOperator>, ExcContractionMap,
                     SymEngine::RCPBasicKeyLess>
        VxcContractionMap;

    // Extract the grid weight, and corresponding contraction between XC
    // functional derivative and generalized density vectors from their
    // multiplication expression.
    std::tuple<SymEngine::RCP<const NonElecFunction>,
               unsigned int,
               SymEngine::RCP<const SymEngine::Basic>>
    extract_exc_contraction(const SymEngine::RCP<const SymEngine::Mul>& expression);

    // Add the output from `extract_exc_contraction()` into an `ExcContractionMap`
    inline void add_exc_contraction(
        ExcContractionMap& energyMap,
        const std::tuple<SymEngine::RCP<const NonElecFunction>,
                         unsigned int,
                         SymEngine::RCP<const SymEngine::Basic>>& value
    )
    {
        // Check if the grid weight exists in the outer map
        auto weight_map = energyMap.find(std::get<0>(value));
        if (weight_map == energyMap.end()) {
            energyMap.emplace(
                std::get<0>(value),
                ExcDensityContractionMap({{std::get<1>(value), std::get<2>(value)}})
            );
        }
        else {
            // Check if the functional derivative of XC energy density exists
            // in the inner map
            auto exc_map = weight_map->second.find(std::get<1>(value));
            if (exc_map == weight_map->second.end()) {
                weight_map->second.emplace(std::get<1>(value), std::get<2>(value));
            }
            else {
                exc_map->second = SymEngine::add(exc_map->second, std::get<2>(value));
            }
        }
    }

    // Extract all terms in XC energy or its derivatives, i.e. (un)perturbed
    // weights, XC functional derivative vectors and perturbed generalized
    // density vectors. Results are arranged in a nested map. The key of the
    // outer map is (un)perturbed weights, and the value of the outer map is
    // another map whose key is the order of functional derivatives of XC
    // energy density, and whose value is the corresponding perturbed
    // generalized density vectors.
    inline ExcContractionMap extract_energy_map(
        const SymEngine::RCP<const SymEngine::Basic>& expression
    )
    {
        // XC energy or its derivatives must be either `SymEngine::Mul` or
        // `SymEngine::Add`
        SYMENGINE_ASSERT(
            SymEngine::is_a<const SymEngine::Mul>(*expression) ||
            SymEngine::is_a<const SymEngine::Add>(*expression)
        )
        // Unperturbed or the first-order case
        if (SymEngine::is_a_sub<const SymEngine::Mul>(*expression)) {
            auto contr_term = extract_exc_contraction(
                SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(expression)
            );
            return ExcContractionMap({
                {
                    std::get<0>(contr_term),
                    ExcDensityContractionMap({
                        {std::get<1>(contr_term), std::get<2>(contr_term)}
                    })
                }
            });
        }
        // Perturbed case
        else {
            ExcContractionMap energy_map;
            auto energy = SymEngine::rcp_dynamic_cast<const SymEngine::Add>(expression);
            auto energy_terms = energy->get_args();
            // No coefficient exists in the XC energy derivatives
            SYMENGINE_ASSERT(
                !SymEngine::is_a_sub<const SymEngine::Number>(*energy_terms.front())
            )
            for (const auto& term: energy_terms) {
                SYMENGINE_ASSERT(
                    SymEngine::is_a_sub<const SymEngine::Mul>(*term)
                )
                auto contr_term = extract_exc_contraction(
                    SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(term)
                );
                add_exc_contraction(energy_map, contr_term);
            }
            return energy_map;
        }
    }

    // Extract grid weight, XC functional derivative, generalized density
    // vectors and overlap distributions from their multiplication expression.
    std::pair<SymEngine::RCP<const OneElecOperator>, ExcContractionMap>
    extract_vxc_contraction(
        const SymEngine::RCP<const SymEngine::MatrixMul>& expression
    );

    // Merge the second `ExcContractionMap` to the first one
    void merge_exc_contraction(ExcContractionMap& map1, const ExcContractionMap& map2);

    // Extract all terms in XC potential operator or its derivatives, i.e.
    // (un)perturbed weights, XC functional derivative vectors, perturbed
    // generalized density vectors and (un)perturbed generalized overlap
    // distributions. Results are arranged in a nested map. The key of the
    // outermost map is (un)perturbed generalized overlap distributions, whose
    // value is `ExcContractionMap`.
    inline VxcContractionMap extract_potential_map(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& expression
    )
    {
        // XC potential or its derivatives must be either
        // `SymEngine::MatrixMul` or `SymEngine::MatrixAdd`
        SYMENGINE_ASSERT(
            SymEngine::is_a<const SymEngine::MatrixMul>(*expression) ||
            SymEngine::is_a<const SymEngine::MatrixAdd>(*expression)
        )
        // Unperturbed case or when the generalized overlap distribution does
        // not depend on the applied perturbation(s)
        if (SymEngine::is_a_sub<const SymEngine::MatrixMul>(*expression)) {
            auto potential_term = extract_vxc_contraction(
                SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(expression)
            );
            return VxcContractionMap({potential_term});
        }
        // Perturbed case
        else {
            VxcContractionMap potential_map;
            auto potential = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixAdd>(expression);
            auto potential_terms = potential->get_args();
            for (const auto& term: potential_terms) {
                SYMENGINE_ASSERT(
                    SymEngine::is_a_sub<const SymEngine::MatrixMul>(*term)
                )
                auto vxc_factors = extract_vxc_contraction(
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(term)
                );
                // Check if the generalized overlap distribution exists
                // in the outermost map
                auto energy_map = potential_map.find(vxc_factors.first);
                if (energy_map == potential_map.end()) {
                    potential_map.emplace(vxc_factors);
                }
                else {
                    merge_exc_contraction(energy_map->second, vxc_factors.second);
                }
            }
            return potential_map;
        }
    }

    // Check if two `ExcContractionMap`'s are equivalent
    bool eq_exc_contraction(
        const ExcContractionMap& map1, const ExcContractionMap& map2
    );

    // Check if two `VxcContractionMap`'s are equivalent
    bool eq_vxc_contraction(
        const VxcContractionMap& map1, const VxcContractionMap& map2
    );
}
