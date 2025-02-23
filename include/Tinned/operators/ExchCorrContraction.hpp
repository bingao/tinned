/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file implements functions for contractions between exchange-correlation
   energy functional derivative vectors and generalized density vectors.

   2024-05-10, Bin Gao:
   * add functions to stringify maps

   2023-11-08, Bin Gao:
   * first version
*/

#pragma once

#include <map>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_rcp.h>

#include "Tinned/operators/NonElecFunction.hpp"
#include "Tinned/operators/ElectronicState.hpp"
#include "Tinned/operators/OneElecOperator.hpp"
#include "Tinned/operators/CompositeFunction.hpp"

#include "Tinned/visitors/StringifyVisitor.hpp"

namespace Tinned
{
    // Type for contractions between XC functional derivative vectors and
    // perturbed generalized density vectors
    typedef std::map<SymEngine::RCP<const CompositeFunction>,
                     SymEngine::RCP<const SymEngine::Basic>,
                     SymEngine::RCPBasicKeyLess>
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
    typedef std::map<SymEngine::RCP<const NonElecFunction>,
                     ExcDensityContractionMap,
                     SymEngine::RCPBasicKeyLess>
        ExcContractionMap;

    // Type for the collection of (un)perturbed weights, XC functional
    // derivative vectors, perturbed generalized density vectors and
    // (un)perturbed generalized overlap distributions
    typedef std::map<SymEngine::RCP<const OneElecOperator>, ExcContractionMap,
                     SymEngine::RCPBasicKeyLess>
        VxcContractionMap;

    // Make (unperturbed) generalized density vector
    //FIXME: change `ElectronicState` to OneElecDensity?
    inline SymEngine::RCP<const SymEngine::Basic> make_density_vector(
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega
    )
    {
        return SymEngine::trace(SymEngine::matrix_mul({Omega, state}));
    }

    // Make XC energy density from electronic state, overlap distribution and
    // the order of XC energy functional
    inline SymEngine::RCP<const CompositeFunction> make_exc_density(
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const unsigned int order
    )
    {
        return SymEngine::make_rcp<const CompositeFunction>(
            // XC energy density as the outer function
            std::string("exc"),
            // Generalized density vector as the inner function
            make_density_vector(state, Omega),
            order
        );
    }

    // Check if two `ExcDensityContractionMap`'s are equivalent
    inline bool eq_exc_density_map(
        const ExcDensityContractionMap& map1, const ExcDensityContractionMap& map2
    )
    {
        if (map1.size()==map2.size()) {
            auto term1 = map1.begin();
            auto term2 = map2.begin();
            for (; term1!=map1.end(); ++term1, ++term2) {
                if (SymEngine::eq(*term1->first, *term2->first)) {
                    if (term1->second.is_null()) {
                        if (!term2->second.is_null()) return false;
                    }
                    else {
                        if (term2->second.is_null()) return false;
                        if (SymEngine::neq(*term1->second, *term2->second)) return false;
                    }
                }
                else {
                    return false;
                }
            }
            return true;
        }
        else {
            return false;
        }
    }

    // Stringify `ExcDensityContractionMap`
    inline std::string stringify_exc_density_map(
        StringifyVisitor& visitor,
        const ExcDensityContractionMap& exc_map,
        const std::string& str_indent = std::string("    ")
    )
    {
        auto str_term_indent = str_indent + std::string("    ");
        std::ostringstream o;
        o << str_indent << "{\n";
        for (const auto& term: exc_map) {
            o << str_term_indent << visitor.convert(*term.first) << ",\n";
            if (term.second.is_null()) {
                o << str_term_indent << "null,\n";
            }
            else {
                o << str_term_indent << visitor.convert(*term.second) << ",\n";
            }
        }
        o << str_indent << "}";
        return o.str();
    }

    // Extract the grid weight, and corresponding contraction between XC
    // functional derivative and generalized density vectors from their
    // multiplication expression.
    std::tuple<SymEngine::RCP<const NonElecFunction>,
               SymEngine::RCP<const CompositeFunction>,
               SymEngine::RCP<const SymEngine::Basic>>
    extract_exc_contraction(const SymEngine::RCP<const SymEngine::Mul>& expr);

    // Add the output from `extract_exc_contraction()` into an `ExcContractionMap`
    void add_exc_contraction(
        ExcContractionMap& energy_map,
        const std::tuple<SymEngine::RCP<const NonElecFunction>,
                         SymEngine::RCP<const CompositeFunction>,
                         SymEngine::RCP<const SymEngine::Basic>>& value
    );

    // Extract all terms in XC energy or its derivatives, i.e. (un)perturbed
    // weights, XC functional derivative vectors and perturbed generalized
    // density vectors. Results are arranged in a nested map. The key of the
    // outer map is (un)perturbed weights, and the value of the outer map is
    // another map whose key is the order of functional derivatives of XC
    // energy density, and whose value is the corresponding perturbed
    // generalized density vectors.
    ExcContractionMap extract_energy_map(
        const SymEngine::RCP<const SymEngine::Basic>& expr
    );

    // Merge the second `ExcContractionMap` to the first one
    void merge_energy_map(ExcContractionMap& map1, const ExcContractionMap& map2);

    // Convert `ExcContractionMap` to XC energy or its derivatives
    SymEngine::RCP<const SymEngine::Basic> convert_energy_map(
        const ExcContractionMap& energy_map
    );

    // Canonicalize the expression of XC energy or its derivatives
    inline SymEngine::RCP<const SymEngine::Basic> canonicalize_xc_energy(
        const SymEngine::RCP<const SymEngine::Basic>& expr
    )
    {
        return convert_energy_map(extract_energy_map(expr));
    }

    // Check if two `ExcContractionMap`'s are equivalent
    bool eq_energy_map(const ExcContractionMap& map1, const ExcContractionMap& map2);

    // Stringify `ExcContractionMap`
    inline std::string stringify_energy_map(
        StringifyVisitor& visitor,
        const ExcContractionMap& energy_map,
        const std::string& str_indent = std::string("    ")
    )
    {
        auto str_term_indent = str_indent + std::string("    ");
        std::ostringstream o;
        o << str_indent << "{\n";
        for (const auto& weight_map: energy_map) {
            o << str_term_indent << visitor.convert(*weight_map.first) << ",\n"
              << stringify_exc_density_map(visitor, weight_map.second, str_term_indent)
              << ",\n";
        }
        o << str_indent << "}";
        return o.str();
    }

    // Extract grid weight, XC functional derivative, generalized density
    // vectors and overlap distributions from their multiplication expression.
    std::pair<SymEngine::RCP<const OneElecOperator>, ExcContractionMap>
    extract_vxc_contraction(
        const SymEngine::RCP<const SymEngine::MatrixMul>& expr
    );

    // Extract all terms in XC potential operator or its derivatives, i.e.
    // (un)perturbed weights, XC functional derivative vectors, perturbed
    // generalized density vectors and (un)perturbed generalized overlap
    // distributions. Results are arranged in a nested map. The key of the
    // outermost map is (un)perturbed generalized overlap distributions, whose
    // value is `ExcContractionMap`.
    VxcContractionMap extract_potential_map(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& expr
    );

    // Merge the second `VxcContractionMap` to the first one
    void merge_potential_map(VxcContractionMap& map1, const VxcContractionMap& map2);

    // Convert `VxcContractionMap` to XC potential operator or its derivatives
    SymEngine::RCP<const SymEngine::MatrixExpr> convert_potential_map(
        const VxcContractionMap& potential_map
    );

    // Canonicalize the expression of XC potential operator or its derivatives
    inline SymEngine::RCP<const SymEngine::MatrixExpr> canonicalize_xc_potential(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& expr
    )
    {
        return convert_potential_map(extract_potential_map(expr));
    }

    // Check if two `VxcContractionMap`'s are equivalent
    bool eq_potential_map(const VxcContractionMap& map1, const VxcContractionMap& map2);

    // Stringify `VxcContractionMap`
    inline std::string stringify_potential_map(
        StringifyVisitor& visitor,
        const VxcContractionMap& potential_map,
        const std::string& str_indent = std::string("    ")
    )
    {
        auto str_term_indent = str_indent + std::string("    ");
        std::ostringstream o;
        o << str_indent << "{\n";
        for (const auto& energy_map: potential_map) {
            o << str_term_indent << visitor.convert(*energy_map.first) << ",\n"
              << stringify_energy_map(visitor, energy_map.second, str_term_indent)
              << ",\n";
        }
        o << str_indent << "}";
        return o.str();
    }
}
