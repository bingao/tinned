/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of finding a given symbol and all its
   differentiated ones.

   2024-09-18, Bin Gao:
   * use `SymEngine::vec_basic` for return result

   2024-05-07, Bin Gao:
   * change `find_all` to a function that returns `SymEngine::set_basic`

   2023-10-28, Bin Gao:
   * first version
*/

#pragma once

#include <functional>
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

#include "Tinned/PertDependency.hpp"

#include "Tinned/PerturbedParameter.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecEnergy.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOverlap.hpp"

namespace Tinned
{
    // Forward declaration
    class ExchCorrEnergy;
    class ExchCorrPotential;

    class FindAllVisitor: public SymEngine::BaseVisitor<FindAllVisitor>
    {
        protected:
            SymEngine::vec_basic result_;
            SymEngine::RCP<const SymEngine::Basic> symbol_;

            // Symbols are sorted according to their derivatives and hash
            template<typename T,
                     typename std::enable_if<
                         std::is_same<T, const SymEngine::MatrixDerivative>::value ||
                         std::is_same<T, const SymEngine::Derivative>::value, int>::type = 0>
            inline void insert_symbol(T& x)
            {
                if (result_.empty()) {
                    result_.push_back(x.rcp_from_this());
                }
                else {
                    auto order_x = x.get_symbols().size();
                    auto hash_x = x.hash();
                    auto iter = result_.begin();
                    for (; iter!=result_.end(); ++iter) {
                        auto s = SymEngine::rcp_dynamic_cast<const T>(*iter);
                        auto order_s = s->get_symbols().size();
                        if (order_x==order_s) {
                            if (SymEngine::eq(x, *s)) return;
                            auto hash_s = s->hash();
                            if (hash_x<hash_s) {
                                result_.insert(iter, x.rcp_from_this());
                                return;
                            }
                        }
                        else if (order_x<order_s) {
                            result_.insert(iter, x.rcp_from_this());
                            return;
                        }
                    }
                    // Add `x` to the end
                    result_.push_back(x.rcp_from_this());
                }
            }

            template<typename T,
                     typename std::enable_if<
                         std::is_same<T, const PerturbedParameter>::value ||
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
                if (result_.empty()) {
                    result_.push_back(x.rcp_from_this());
                }
                else {
                    auto order_x = x.get_derivatives().size();
                    auto hash_x = x.hash();
                    auto iter = result_.begin();
                    for (; iter!=result_.end(); ++iter) {
                        auto s = SymEngine::rcp_dynamic_cast<const T>(*iter);
                        auto order_s = s->get_derivatives().size();
                        if (order_x==order_s) {
                            if (SymEngine::eq(x, *s)) return;
                            auto hash_s = s->hash();
                            if (hash_x<hash_s) {
                                result_.insert(iter, x.rcp_from_this());
                                return;
                            }
                        }
                        else if (order_x<order_s) {
                            result_.insert(iter, x.rcp_from_this());
                            return;
                        }
                    }
                    result_.push_back(x.rcp_from_this());
                }
            }

            template<typename T,
                     typename std::enable_if<
                         !(std::is_same<T, const SymEngine::MatrixDerivative>::value ||
                         std::is_same<T, const SymEngine::Derivative>::value ||
                         std::is_same<T, const PerturbedParameter>::value ||
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
                if (result_.empty()) {
                    result_.push_back(x.rcp_from_this());
                }
                else {
                    auto hash_x = x.hash();
                    auto iter = result_.begin();
                    for (; iter!=result_.end(); ++iter) {
                        auto hash_s = (*iter)->hash();
                        if (hash_x==hash_s) {
                            if (SymEngine::eq(x, *(*iter))) return;
                        }
                        else if (hash_x<hash_s) {
                            result_.insert(iter, x.rcp_from_this());
                            return;
                        }
                    }
                    result_.push_back(x.rcp_from_this());
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

            // Method called by objects to prcess their argument(s)
            inline void apply_(const SymEngine::RCP<const SymEngine::Basic>& x)
            {
                x->accept(*this);
            }

        public:
            explicit FindAllVisitor(
                const SymEngine::RCP<const SymEngine::Basic>& symbol
            ) : symbol_(symbol) {}

            inline SymEngine::vec_basic apply(
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
    inline SymEngine::vec_basic find_all(
        const SymEngine::RCP<const SymEngine::Basic>& x,
        const SymEngine::RCP<const SymEngine::Basic>& symbol
    )
    {
        FindAllVisitor visitor(symbol);
        return visitor.apply(x);
    }
}
