/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of finding a given symbol and all its
   differentiated ones.

   2024-11-10, Bin Gao:
   * use `FindAllResult` for result which is well organized according to the
     order of derivatives

   2024-09-18, Bin Gao:
   * use `SymEngine::vec_basic` for return result

   2024-05-07, Bin Gao:
   * change `find_all` to a function that returns `SymEngine::set_basic`

   2023-10-28, Bin Gao:
   * first version
*/

#pragma once

#include <functional>
#include <map>
#include <type_traits>

#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/mul.h>
#include <symengine/number.h>
#include <symengine/symbol.h>
#include <symengine/matrices/conjugate_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_derivative.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/trace.h>
#include <symengine/matrices/transpose.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/symengine_rcp.h>
#include <symengine/visitor.h>

#include "Tinned/operators/PerturbedParameter.hpp"
#include "Tinned/operators/ElectronicState.hpp"
#include "Tinned/operators/OneElecDensity.hpp"
#include "Tinned/operators/OneElecOperator.hpp"
#include "Tinned/operators/TwoElecEnergy.hpp"
#include "Tinned/operators/TwoElecOperator.hpp"
#include "Tinned/operators/NonElecFunction.hpp"
#include "Tinned/operators/TemporumOverlap.hpp"

namespace Tinned
{
    // Forward declaration
    class ExchCorrEnergy;
    class ExchCorrPotential;

    // Type for result from `FindAllVisitor`
    typedef std::map<unsigned int, SymEngine::set_basic> FindAllResult;

    class FindAllVisitor: public SymEngine::BaseVisitor<FindAllVisitor>
    {
        protected:
            FindAllResult result_;
            SymEngine::RCP<const SymEngine::Basic> symbol_;

            // Symbols are sorted according to their derivatives and hash
            template<typename T,
                     typename std::enable_if<
                         std::is_same<T, const SymEngine::MatrixDerivative>::value ||
                         std::is_same<T, const SymEngine::Derivative>::value, int>::type = 0>
            inline void insert_symbol(T& x)
            {
                auto order = x.get_symbols().size();
                auto iter = result_.find(order);
                if (iter==result_.end()) {
                    result_.insert({order, SymEngine::set_basic({x.rcp_from_this()})});
                }
                else {
                    iter->second.insert(x.rcp_from_this());
                }
            }

            inline void insert_symbol(const CompositeFunction& x)
            {
                auto order = x.get_order();
                auto iter = result_.find(order);
                if (iter==result_.end()) {
                    result_.insert({order, SymEngine::set_basic({x.rcp_from_this()})});
                }
                else {
                    iter->second.insert(x.rcp_from_this());
                }
            }

            template<typename T,
                     typename std::enable_if<
                         std::is_same<T, const PerturbedParameter>::value ||
                         std::is_same<T, const ElectronicState>::value ||
                         std::is_same<T, const OneElecDensity>::value ||
                         std::is_same<T, const OneElecOperator>::value ||
                         std::is_same<T, const TwoElecEnergy>::value ||
                         std::is_same<T, const TwoElecOperator>::value ||
                         std::is_same<T, const ExchCorrEnergy>::value ||
                         std::is_same<T, const ExchCorrPotential>::value ||
                         std::is_same<T, const NonElecFunction>::value ||
                         std::is_same<T, const TemporumOverlap>::value, int>::type = 0>
            inline void insert_symbol(T& x)
            {
                auto order = x.get_derivatives().size();
                auto iter = result_.find(order);
                if (iter==result_.end()) {
                    result_.insert({order, SymEngine::set_basic({x.rcp_from_this()})});
                }
                else {
                    iter->second.insert(x.rcp_from_this());
                }
            }

            template<typename T,
                     typename std::enable_if<
                         !(std::is_same<T, const SymEngine::MatrixDerivative>::value ||
                         std::is_same<T, const SymEngine::Derivative>::value ||
                         std::is_same<T, const CompositeFunction>::value ||
                         std::is_same<T, const PerturbedParameter>::value ||
                         std::is_same<T, const ElectronicState>::value ||
                         std::is_same<T, const OneElecDensity>::value ||
                         std::is_same<T, const OneElecOperator>::value ||
                         std::is_same<T, const TwoElecEnergy>::value ||
                         std::is_same<T, const TwoElecOperator>::value ||
                         std::is_same<T, const ExchCorrEnergy>::value ||
                         std::is_same<T, const ExchCorrPotential>::value ||
                         std::is_same<T, const NonElecFunction>::value ||
                         std::is_same<T, const TemporumOverlap>::value), int>::type = 0>
            inline void insert_symbol(T& x)
            {
                auto iter = result_.find(0);
                if (iter==result_.end()) {
                    result_.insert({0, SymEngine::set_basic({x.rcp_from_this()})});
                }
                else {
                    iter->second.insert(x.rcp_from_this());
                }
            }

            // Function template to check if `x` is the symbol we want to find
            // according to a given comparision condition, and update `result_`
            // when `x` is the symbol to be found
            template<typename T> inline bool find_with_condition(
                T& x,
                const std::function<bool(T&, T&)>& condition
            )
            {
                if (SymEngine::is_a_sub<T>(*symbol_)) {
                    auto& s = SymEngine::down_cast<T&>(*symbol_);
                    if (condition(x, s)) {
                        insert_symbol(x);
                        return true;
                    }
                }
                return false;
            }

            // Function template for objects that requires equivalence comparison
            template<typename T> inline bool find_equivalence(T& x)
            {
                if (symbol_->__eq__(x)) {
                    insert_symbol(x);
                    return true;
                }
                return false;
            }

            // Function template for objects that compares only the names
            template<typename T> inline bool find_only_name(T& x)
            {
                return find_with_condition<T>(
                    x,
                    [&](T& op1, T& op2) -> bool { return op1.get_name()==op2.get_name(); }
                );
            }

            // Function template for objects without arguments, we compare only
            // their names and dependencies
            template<typename T> inline void find_with_dependencies(T& x)
            {
                find_with_condition<T>(
                    x,
                    [&](T& op1, T& op2) -> bool {
                        return op1.get_name()==op2.get_name() &&
                            eq_dependency(op1.get_dependencies(), op2.get_dependencies());
                    }
                );
            }

            // Function template for one argument function like classes and the
            // name of the function does not matter
            template<typename Fun, typename Arg>
            inline void find_one_arg_f(
                Fun& x,
                const std::function<bool(Fun&, Fun&)>& condition,
                const SymEngine::RCP<Arg>& arg
            )
            {
                if (!find_with_condition<Fun>(x, condition)) apply_(arg);
            }

            // Function to check the name and dependencies of `TwoElecOperator`
            // objects and the name of their denisty matrices
            inline bool comp_2el_operator(
                const TwoElecOperator& op1, const TwoElecOperator& op2
            ) const
            {
                return op1.get_name()==op2.get_name()
                    && eq_dependency(op1.get_dependencies(), op2.get_dependencies())
                    && op1.get_state()->get_name()==op2.get_state()->get_name();
            }

            // Method called by objects to process their argument(s)
            inline void apply_(const SymEngine::RCP<const SymEngine::Basic>& x)
            {
                x->accept(*this);
            }

        public:
            explicit FindAllVisitor(
                const SymEngine::RCP<const SymEngine::Basic>& symbol
            ) : symbol_(symbol) {}

            inline FindAllResult apply(
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

    // Helper function to find a given `symbol` and all its differentiated ones in `x`
    inline FindAllResult find_all(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::RCP<const SymEngine::Basic>& symbol
    )
    {
        FindAllVisitor visitor(symbol);
        return visitor.apply(x);
    }

    // Equality comparison of two found results
    inline bool is_equal(const FindAllResult& lhs, const FindAllResult& rhs)
    {
        if (lhs.size()!=rhs.size()) return false;
        auto iter_lhs = lhs.begin();
        auto iter_rhs = rhs.begin();
        for (; iter_lhs!=lhs.end(); ++iter_lhs,++iter_rhs) {
            if (iter_lhs->first!=iter_rhs->first) return false;
            if (!SymEngine::unified_eq(iter_lhs->second, iter_rhs->second)) return false;
        }
        return true;
    }
}
