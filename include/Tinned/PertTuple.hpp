/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of perturbation tuples.

   2024-04-30, Bin Gao:
   * first version
*/

#pragma once

#include <iostream>
#include "Tinned/StringifyVisitor.hpp"

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/symbol.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_rcp.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/ZeroOperator.hpp"
#include "Tinned/ZerosRemover.hpp"

namespace Tinned
{
    // Type for perturbation tuple -- an ordered list of the perturbation strengths
    typedef std::multiset<SymEngine::RCP<const Perturbation>, SymEngine::RCPBasicKeyLess>
        PertTuple;

    // Convert derivatives to `PertTuple`
    inline PertTuple make_pert_tuple(const SymEngine::multiset_basic& derivatives)
    {
        PertTuple perturbations;
        for (const auto& p: derivatives) {
            SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(*p))
            perturbations.insert(SymEngine::rcp_dynamic_cast<const Perturbation>(p));
        }
        return perturbations;
    }

    // Helper function to do high-order differentiation, and to remove zero
    // quantities
    template<typename T,
             typename std::enable_if<std::is_same<T, PertTuple>::value ||
                 std::is_same<T, SymEngine::multiset_basic>::value ||
                 std::is_same<T, SymEngine::vec_basic>::value, int>::type = 0>
    inline SymEngine::RCP<const SymEngine::Basic> differentiate(
        const SymEngine::RCP<const SymEngine::Basic>& expr, const T& perturbations
    )
    {
        auto result = expr;
        for (const auto& p: perturbations) {
            SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(*p))
            result = result->diff(SymEngine::rcp_dynamic_cast<const Perturbation>(p));
        }
        result = remove_zeros(result);
        if (result.is_null()) {
            if (SymEngine::is_a_sub<const SymEngine::MatrixExpr>(*expr)) {
                return make_zero_operator();
            }
            else {
                return SymEngine::zero;
            }
        }
        else {
            return result;
        }
    }

    // Compute the sum of perturbation frequencies
    template<typename T,
             typename std::enable_if<std::is_same<T, PertTuple>::value ||
                 std::is_same<T, SymEngine::multiset_basic>::value ||
                 std::is_same<T, SymEngine::vec_basic>::value, int>::type = 0>
    inline SymEngine::RCP<const SymEngine::Basic> get_frequency_sum(
        const T& perturbations
    )
    {
        SymEngine::RCP<const SymEngine::Basic> result = SymEngine::zero;
        for (const auto& p: perturbations) {
            SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(*p))
            result = SymEngine::add(
                result,
                SymEngine::rcp_dynamic_cast<const Perturbation>(p)->get_frequency()
            );
        }
        return result;
    }

    // For a set of perturbations {b1, b2, ..., bn} and an `order` of
    // perturbation-strength derivatives, we first generate a perturbation
    // tuple of the `order`, and then build all its perturbation-strength
    // derivatives by permuting its perturbations.
    //
    // Examples of the 3rd order perturbation tuples from different sets
    // {a}, 3rd order: aaa
    // {a, b}, 3rd order: aaa (3+0), aab (2+1), abb (1+2), bbb (0+3)
    // {a, b, c}, 3rd order: aaa (3+0+0), aab (2+1+0), aac (2+0+1),
    //                       abb (1+2+0), abc (1+1+1), acc (1+0+2),
    //                       bbb (0+3+0), bbc (0+2+1), bcc (0+1+2), ccc (0+0+3)
    //
    // So, the generation of perturbation tuples is equivalent to find all weak
    // compositions of the `order`. The number of compositions is equal to the
    // size of the set of perturbations. It is also the method "Functions from
    // N to X, up to a permutation of N" in Twelvefold way. (2024, February 4),
    // in Wikipedia, https://en.wikipedia.org/wiki/Twelvefold_way
    //
    // The algorithm NEXCOM, in Albert Nijenhuis and Herbert S.  Wilf,
    // "Combinatorial Algorithms: For Computers and Calculators", 2nd edition
    // (1978) can be used for the generation of perturbation tuples. But we
    // modify that from https://stackoverflow.com/a/17463867 in the current
    // implementation.
    class PertPermutation
    {
        protected:
            unsigned int order_;
            SymEngine::set_basic perturbations_;
            std::vector<unsigned int> pert_compositions_;
            std::size_t iter_compositions_;
            // Indicate if initializing perturbation tuple for the first time use
            bool init_pert_tuple_;
            SymEngine::vec_basic pert_tuple_;
            bool pert_tuple_remaining_;
            std::vector<unsigned int> permut_positions_;

            // Move to next weak composition, modified from https://stackoverflow.com/a/17463867
            inline bool next_weak_compositions() noexcept
            {
                if (iter_compositions_+1!=pert_compositions_.size()) {
                    --pert_compositions_[iter_compositions_];
                    ++iter_compositions_;
                    pert_compositions_[iter_compositions_] = 1;
                }
                else {
                    // Reach the last weak composition and reset to the first one
                    if (pert_compositions_.back()==order_) {
                        std::fill(pert_compositions_.begin(), pert_compositions_.end(), 0);
                        pert_compositions_[0] = order_;
                        iter_compositions_ = 0;
                        return false;
                    }
                    else {
                        --iter_compositions_;
                        for (; pert_compositions_[iter_compositions_]==0; --iter_compositions_);
                        --pert_compositions_[iter_compositions_];
                        unsigned int tmp = pert_compositions_.back()+1;
                        ++iter_compositions_;
                        pert_compositions_.back() = 0;
                        pert_compositions_[iter_compositions_] = tmp;
                    }
                }
                return true;
            }

        public:
            explicit PertPermutation(
                const unsigned int order,
                const SymEngine::set_basic& perturbations
            ) noexcept: order_(order), perturbations_(perturbations)
            {
                pert_compositions_ = std::vector<unsigned int>(perturbations.size(), 0);
                pert_compositions_[0] = order_;
                iter_compositions_ = 0;
                init_pert_tuple_ = true;
                permut_positions_.reserve(order);
                for (unsigned int i=0; i<order; ++i) permut_positions_.push_back(i);
            }

            // Get the next perturbation tuple
            inline SymEngine::vec_basic get_pert_tuple(bool& remaining) noexcept
            {
                auto pert_tuple = SymEngine::vec_basic({});
                auto iter_pert = perturbations_.begin();
                for (std::size_t p=0; p<pert_compositions_.size(); ++p,++iter_pert)
                    if (pert_compositions_[p]>0)
                        for (std::size_t i=0; i<pert_compositions_[p]; ++i)
                            pert_tuple.push_back(*iter_pert);
                if (perturbations_.size()==1) {
                    remaining = false;
                }
                else {
                    remaining = next_weak_compositions();
                }
                return pert_tuple;
            }

            // Get the next permuting perturbation-strength derivatives
            inline SymEngine::vec_basic get_derivatives(bool& remaining) noexcept
            {
                if (init_pert_tuple_) {
                    pert_tuple_ = get_pert_tuple(pert_tuple_remaining_);
                    init_pert_tuple_ = false;
                }
                SymEngine::vec_basic derivatives;
                for (const auto& p: permut_positions_) derivatives.push_back(pert_tuple_[p]);
                // Move to the next permuting perturbation-strength derivative
                remaining = std::next_permutation(
                    permut_positions_.begin(), permut_positions_.end()
                );
                if (!remaining) {
                    if (pert_tuple_remaining_) remaining = true;
                    pert_tuple_ = get_pert_tuple(pert_tuple_remaining_);
                }
                return derivatives;
            }

            ~PertPermutation() noexcept = default;
    };
}
