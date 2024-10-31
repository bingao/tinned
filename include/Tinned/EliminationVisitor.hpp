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
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/number.h>
#include <symengine/mul.h>
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
#include "Tinned/PertMultichain.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/PerturbedParameter.hpp"

#include "Tinned/VisitorUtilities.hpp"

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
                         std::is_same<T, const PerturbedParameter>::value, int>::type = 0>
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

            // Function template to eliminate a (response) parameter
            template<typename T> inline void eliminate_parameter(T& x)
            {
                result_ = is_parameter_eliminable(x)
                    ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
            }

            // Function template for a function like object with one or more arguments
            template<typename Fun, typename FirstArg, typename... Args>
            inline void eliminate_a_function(
                const Fun& x,
                const std::function<SymEngine::RCP<const SymEngine::Basic>(
                    const SymEngine::vec_basic&
                )>& constructor,
                const FirstArg& first_arg,
                const Args&... args
            )
            {
                // We check only if each argument has parameter(s) to be eliminated
                auto f_args = SymEngine::vec_basic({});
                auto has_arg_affected = false;
                // `result_` will be null if any argument is eliminated completely
                if (visit_arguments(
                    *this, f_args, has_arg_affected, first_arg, args...
                )) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    result_ = has_arg_affected ? constructor(f_args) : x.rcp_from_this();
                }
            }

        public:
            explicit EliminationVisitor(
                const SymEngine::RCP<const SymEngine::Basic>& parameter,
                const PertMultichain& perturbations,
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
            void bvisit(const SymEngine::Number& x);
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
    // half length of `perturbations`, and for multipliers, it should be greater
    // than or equal to the ceiling function of the half length of
    // `perturbations` according to J. Chem. Phys. 129, 214103 (2008).
    inline SymEngine::RCP<const SymEngine::Basic> eliminate(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::RCP<const SymEngine::Basic>& parameter,
        const PertMultichain& perturbations,
        const unsigned int min_order
    )
    {
        EliminationVisitor visitor(parameter, perturbations, min_order);
        return visitor.apply(x);
    }
}
