/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of elimination of response parameters by
   following J. Chem. Phys. 129, 214103 (2008).

   2024-05-03, Bin Gao:
   * previous implementation does not work since it uses `SymEngine::eq()` to
     compare a symbol and the parameter to be eliminated, which will also
     compare their derivatives.

   2024-04-26, Bin Gao:
   * first version
*/

#pragma once

#include <functional>
#include <string>
#include <type_traits>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/complex.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/rational.h>
#include <symengine/symbol.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/StateVector.hpp"
#include "Tinned/LagMultiplier.hpp"

namespace Tinned
{
    class EliminationVisitor: public SymEngine::BaseVisitor<EliminationVisitor>
    {
        protected:
            SymEngine::RCP<const SymEngine::Basic> result_;
            SymEngine::RCP<const SymEngine::Basic> parameter_;
            SymEngine::set_basic perturbations_;
            unsigned int max_order_;
            unsigned int min_order_;

            // Check the order of derivatives
            inline bool match_derivatives(const SymEngine::multiset_basic& derivatives) const
            {
                unsigned int order = 0;
                for (const auto& p: perturbations_) order += derivatives.count(p);
                return (order<=max_order_ && order>=min_order_) ? true : false;
            }

            // Check if a (response) parameter `x` should be eliminated
            template<typename T,
                     typename std::enable_if<
                         std::is_same<T, const OneElecDensity>::value ||
                         std::is_same<T, const StateVector>::value ||
                         std::is_same<T, const LagMultiplier>::value, int>::type = 0>
            inline bool is_parameter_eliminable(T& x) const
            {
                if (x.is_same_parameter(parameter_)) {
                    return match_derivatives(x.get_derivatives());
                }
                else {
                    return false;
                }
            }

            // Check if a wave function parameter `x` should be eliminated
            inline bool is_parameter_eliminable(
                const SymEngine::RCP<const ElectronicState>& x
            ) const
            {
                if (x->is_same_parameter(parameter_)) {
                    return match_derivatives(x->get_derivatives());
                }
                else {
                    return false;
                }
            }

            // Template method to eliminate a (response) parameter
            template<typename T> inline void eliminate_parameter(T& x)
            {
                result_ = is_parameter_eliminable(x)
                    ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
            }

            // Template method for one argument function like classes
            template<typename Fun, typename Arg>
            inline void eliminate_one_arg_f(
                Fun& x,
                const SymEngine::RCP<Arg>& arg,
                std::function<SymEngine::RCP<const SymEngine::Basic>(
                    const SymEngine::RCP<Arg>&
                )> constructor
            )
            {
                // We check only if the argument will be eliminated
                auto new_arg = apply(arg);
                if (new_arg.is_null()) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    if (SymEngine::eq(*arg, *new_arg)) {
                        result_ = x.rcp_from_this();
                    }
                    else {
                        result_ = constructor(SymEngine::rcp_dynamic_cast<Arg>(new_arg));
                    }
                }
            }

        public:
            explicit EliminationVisitor(
                const SymEngine::RCP<const SymEngine::Basic>& parameter,
                const PerturbationTuple& perturbations,
                const unsigned int min_order
            ) : parameter_(parameter),
                perturbations_(SymEngine::set_basic(perturbations.begin(), perturbations.end())),
                max_order_(perturbations.size()),
                min_order_(min_order) {}

            inline SymEngine::RCP<const SymEngine::Basic> apply(
                const SymEngine::RCP<const SymEngine::Basic>& x
            )
            {
                x->accept(*this);
                return result_;
            }

            void bvisit(const SymEngine::Basic& x);
            void bvisit(const SymEngine::Symbol& x);
            void bvisit(const SymEngine::Integer& x);
            void bvisit(const SymEngine::Rational& x);
            void bvisit(const SymEngine::Complex& x);
            void bvisit(const SymEngine::Add& x);
            void bvisit(const SymEngine::Mul& x);
            void bvisit(const SymEngine::Constant& x);
            void bvisit(const SymEngine::FunctionSymbol& x);
            void bvisit(const SymEngine::ZeroMatrix& x);
            void bvisit(const SymEngine::MatrixSymbol& x);
            void bvisit(const SymEngine::Trace& x);
            void bvisit(const SymEngine::ConjugateMatrix& x);
            void bvisit(const SymEngine::Transpose& x);
            void bvisit(const SymEngine::MatrixAdd& x);
            void bvisit(const SymEngine::MatrixMul& x);
            void bvisit(const SymEngine::MatrixDerivative& x);
    };

    // Helper function to eliminate a given response `parameter`'s derivatives
    // from `x`. Maximum order of derivatives to be eliminated is the length of
    // `perturbations`, and minimum order is specified by `min_order`. For wave
    // function parameters, it should be greater than the floor function of the
    // half length of `perturbations`, and or multipliers, it should be greater
    // than or equal to the ceiling function of the half length of
    // `perturbations` according to J. Chem. Phys. 129, 214103 (2008).
    inline SymEngine::RCP<const SymEngine::Basic> eliminate(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::RCP<const SymEngine::Basic>& parameter,
        const PerturbationTuple& perturbations,
        const unsigned int min_order
    )
    {
        EliminationVisitor visitor(parameter, perturbations, min_order);
        return visitor.apply(x);
    }
}
