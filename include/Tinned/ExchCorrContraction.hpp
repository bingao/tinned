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
#include <symengine/mul.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/symengine_rcp.h>

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

    // Extract grid weight, XC functional derivative, generalized density
    // vectors and overlap distributions from their multiplication expression.
    std::pair<SymEngine::RCP<const OneElecOperator>, ExcContractionMap>
    extract_vxc_contraction(
        const SymEngine::RCP<const SymEngine::MatrixMul>& expression
    );

    // Add the output from `extract_exc_contraction()` into an `ExcContractionMap`
    inline void add_exc_contraction(
        ExcContractionMap& contrMap,
        const std::tuple<SymEngine::RCP<const NonElecFunction>,
                         unsigned int,
                         SymEngine::RCP<const SymEngine::Basic>>& value
    )
    {
        // Check if the grid weight exists in the outer map
        auto term = contrMap.find(std::get<0>(value));
        if (term == contrMap.end()) {
            contrMap.emplace(
                std::get<0>(value),
                ExcDensityContractionMap({{std::get<1>(value), std::get<2>(value)}})
            );
        }
        else {
            // Check if the functional derivative of XC energy density exists
            // in the inner map
            auto contr = term->second.find(std::get<1>(value));
            if (contr == term->second.end()) {
                term->second.emplace(std::get<1>(value), std::get<2>(value));
            }
            else {
                contr->second = SymEngine::add(contr->second, std::get<2>(value));
            }
        }
    }

    // Merge the second `ExcContractionMap` to the first one
    void merge_exc_contraction(ExcContractionMap& map1, const ExcContractionMap& map2);

    // Check if two `ExcContractionMap`'s are equivalent
    bool eq_exc_contraction(
        const ExcContractionMap& map1, const ExcContractionMap& map2
    );
}
