/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of utilities for different visitors.

   2024-06-05, Bin Gao:
   * first version
*/

#pragma once

#include <string>
#include <type_traits>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/transpose.h>

#include "Tinned/perturbations/PertDependency.hpp"

#include "Tinned/operators/ConjugateTranspose.hpp"
#include "Tinned/operators/ElectronicState.hpp"
#include "Tinned/operators/CompositeFunction.hpp"
#include "Tinned/operators/NonElecFunction.hpp"
#include "Tinned/operators/TwoElecEnergy.hpp"
#include "Tinned/operators/ExchCorrEnergy.hpp"
#include "Tinned/operators/OneElecOperator.hpp"
#include "Tinned/operators/TwoElecOperator.hpp"
#include "Tinned/operators/ExchCorrPotential.hpp"
#include "Tinned/operators/TemporumOperator.hpp"
#include "Tinned/operators/AdjointMap.hpp"
#include "Tinned/operators/ClusterConjHamiltonian.hpp"

namespace Tinned
{
    // Function template to visit only one argument.
    template<typename Visitor, typename Arg,
             typename std::enable_if<!std::is_same<Arg, SymEngine::vec_basic>::value, int>::type = 0>
    inline bool visit_arguments(
        Visitor& v,
        SymEngine::vec_basic& f_args,
        bool& has_arg_affected,
        const Arg& arg
    )
    {
        auto new_arg = v.apply(arg);
        if (new_arg.is_null()) return true;
        f_args.push_back(new_arg);
        if (SymEngine::neq(*arg, *new_arg)) has_arg_affected = true;
        return false;
    }

    // Function template to visit only one argument of type `SymEngine::vec_basic`.
    template<typename Visitor, typename Arg,
             typename std::enable_if<std::is_same<Arg, SymEngine::vec_basic>::value, int>::type = 0>
    inline bool visit_arguments(
        Visitor& v,
        SymEngine::vec_basic& f_args,
        bool& has_arg_affected,
        const Arg& arg
    )
    {
        for (const auto& term: arg)
            if (visit_arguments(v, f_args, has_arg_affected, term)) return true;
        return false;
    }

    // Function template to visit one or more arguments. `f_args` holds all
    // arguments, either affected or unaffected after visiting.
    // `has_arg_affected` indicates if one or more arguments are affected due
    // to the visiting.
    template<typename Visitor, typename FirstArg, typename... Args>
    inline bool visit_arguments(
        Visitor& v,
        SymEngine::vec_basic& f_args,
        bool& has_arg_affected,
        const FirstArg& first_arg, const Args&... args
    )
    {
        if (visit_arguments(v, f_args, has_arg_affected, first_arg)) return true;
        if (visit_arguments(v, f_args, has_arg_affected, args...)) return true;
        return false;
    }

    // Various inline operator constructors, mostly used by different visitors

    inline SymEngine::RCP<const SymEngine::Basic> construct_trace(
        const SymEngine::vec_basic& args
    )
    {
        return SymEngine::trace(
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0])
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_conjugate_matrix(
        const SymEngine::vec_basic& args
    )
    {
        return SymEngine::conjugate_matrix(
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0])
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_transpose(
        const SymEngine::vec_basic& args
    )
    {
        return SymEngine::transpose(
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0])
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_conjugate_transpose(
        const SymEngine::vec_basic& args
    )
    {
        ConjugateTransposeVisitor visitor;
        return visitor.apply(
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0])
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_composite_function(
        const SymEngine::vec_basic& args,
        const std::string& name,
        const unsigned int order
    )
    {
        return SymEngine::make_rcp<const CompositeFunction>(name, args[0], order);
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_2el_energy(
        const SymEngine::vec_basic& args
    )
    {
        return SymEngine::make_rcp<const TwoElecEnergy>(
            SymEngine::rcp_dynamic_cast<const TwoElecOperator>(args[0]),
            SymEngine::rcp_dynamic_cast<const ElectronicState>(args[1])
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_xc_energy(
        const SymEngine::vec_basic& args,
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const NonElecFunction>& weight
    )
    {
        return SymEngine::make_rcp<const ExchCorrEnergy>(
            name, state, Omega, weight, args[0]
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_2el_operator(
        const SymEngine::vec_basic& args,
        const std::string& name,
        const PertDependency& dependencies,
        const SymEngine::multiset_basic& derivatives
    )
    {
        return SymEngine::make_rcp<const TwoElecOperator>(
            name,
            SymEngine::rcp_dynamic_cast<const ElectronicState>(args[0]),
            dependencies,
            derivatives
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_xc_potential(
        const SymEngine::vec_basic& args,
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const NonElecFunction>& weight
    )
    {
        return SymEngine::make_rcp<const ExchCorrPotential>(
            name,
            state,
            Omega,
            weight,
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0])
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_dt_operator(
        const SymEngine::vec_basic& args,
        const TemporumType type
    )
    {
        return SymEngine::make_rcp<const TemporumOperator>(
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0]),
            type
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_adjoint_map(
        const SymEngine::vec_basic& args
    )
    {
        return SymEngine::make_rcp<const AdjointMap>(
            SymEngine::vec_basic(args.begin(), args.end()-1), args.back()
        );
    }

    inline SymEngine::RCP<const SymEngine::Basic> construct_cc_hamiltonian(
        const SymEngine::vec_basic& args
    )
    {
        return SymEngine::make_rcp<const ClusterConjHamiltonian>(
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0]), args[1]
        );
    }
}
