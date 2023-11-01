/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file contains common definitions and functions for derivatives with
   respect to perturbations.

   2023-09-18, Bin Gao:
   * first version
*/

#pragma once

#include <set>
#include <utility>

#include <symengine/basic.h>
#include <symengine/integer.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/zero_matrix.h>

#include "Tinned/Perturbation.hpp"

namespace Tinned
{
    // Comparison function for std::pair<SymEngine::RCP<const Perturbation>,unsigned int>
    struct PertDependencyKeyLess {
        //! true if `x < y`, false otherwise
        bool operator()(
            const std::pair<SymEngine::RCP<const Perturbation>,unsigned int>& x,
            const std::pair<SymEngine::RCP<const Perturbation>,unsigned int>& y
        ) const
        {
            auto result = SymEngine::unified_compare(x, y);
            if (result == 0) return false;
            return result == -1;
        }
    };

    // Type for perturbations that an operator depends on and their maximum
    // orders that can be differentiated
    typedef std::set<std::pair<SymEngine::RCP<const Perturbation>,unsigned int>,
                     PertDependencyKeyLess>
        PertDependency;

    // Equality comparator for perturbation dependencies
    inline bool eq_dependency(const PertDependency& dep1, const PertDependency& dep2)
    {
        if (dep1.size() != dep2.size()) return false;
        auto p1 = dep1.begin();
        auto p2 = dep2.begin();
        for (; p1 != dep1.end(); ++p1, ++p2) {
            if (not (
                SymEngine::unified_eq(p1->first, p2->first) and
                SymEngine::unified_eq(p1->second, p2->second)
            )) return false;
        }
        return true;
    }

    // Find a perturbation in the dependencies and return its maximum order of
    // differentiation
    inline unsigned int find_dependency(
        const PertDependency& dependencies,
        const SymEngine::RCP<const SymEngine::Symbol>& s
    )
    {
        if (SymEngine::is_a<const Perturbation>(*s)) {
            auto p = SymEngine::rcp_static_cast<const Perturbation>(s);
            for (auto& dep: dependencies) {
                if (dep.first->__eq__(*p)) return dep.second;
            }
            return 0;
        }
        else {
            return 0;
        }
    }

    // Convert the dependencies into vec_basic
    // Should be tested for use
    inline SymEngine::vec_basic to_vec_basic(const PertDependency& dependencies)
    {
        SymEngine::vec_basic args;
        for (auto& dep: dependencies) {
            args.push_back(dep.first);
            args.push_back(SymEngine::integer(dep.second));
        }
        return args;
    }

    // Helper function to return a 1x1 zero matrix, can be used for zero
    // derivative
    inline SymEngine::RCP<const SymEngine::MatrixExpr> make_zero_operator()
    {
        return SymEngine::zero_matrix(SymEngine::integer(1), SymEngine::integer(1));
    }
}
