/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of T matrix.

   2023-10-28, Bin Gao:
   * remove member method get_args()

   2023-10-09, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/number.h>
#include <symengine/integer.h>
#include <symengine/symbol.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/matrix_symbol.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/TemporumOperator.hpp"

namespace Tinned
{
    // T = -1/2(<i\frac{\partial}{\partial t}\chi_{\kappa}|\chi_{\lambda}>
    //   + <\chi_{\kappa}|i\frac{\partial}{\partial t}\chi_{\lambda}>)
    class TemporumOverlap: public SymEngine::MatrixSymbol
    {
        protected:
            // Sum of half time-differentiated bra and ket products
            SymEngine::RCP<const SymEngine::Basic> braket_;

        public:
            explicit TemporumOverlap(const PertDependency& dependencies);
            // Only used for `diff_impl()`
            explicit TemporumOverlap(
                const SymEngine::RCP<const SymEngine::Basic>& braket
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;
            //SymEngine::vec_basic get_args() const override;

            // Override the defaut behaviour for diff
            SymEngine::RCP<const SymEngine::Basic> diff_impl(
                const SymEngine::RCP<const SymEngine::Symbol>& s
            ) const override;

            // Get dependencies
            inline PertDependency get_dependencies() const
            {
                auto target = SymEngine::rcp_dynamic_cast<const NonElecFunction>(
                    std::get<1>(get_braket_product(0))->get_target()
                );
                return target->get_dependencies();
            }

            // Get number of half time-differentiated bra and ket products
            inline std::size_t size() const
            {
                if (SymEngine::is_a_sub<const SymEngine::MatrixAdd>(*braket_)) {
                    auto braket = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixAdd>(braket_);
                    return braket->get_args().size();
                }
                else {
                    return 1;
                }
            }

            // Get sum of half time-differentiated bra and ket products
            inline SymEngine::RCP<const SymEngine::Basic> get_braket() const
            {
                return braket_;
            }

            // Get a half time-differentiated bra and ket product
            inline std::tuple<SymEngine::RCP<const SymEngine::Number>,
                              SymEngine::RCP<const TemporumOperator>,
                              SymEngine::RCP<const TemporumOperator>>
            get_braket_product(const std::size_t index) const
            {
                SYMENGINE_ASSERT(index < size())
                auto product = SymEngine::is_a_sub<const SymEngine::MatrixAdd>(*braket_)
                    ? SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>((
                          SymEngine::rcp_dynamic_cast<const SymEngine::MatrixAdd>(
                              braket_
                          )->get_args()
                      )[index])
                    : SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(braket_);

                auto coef = SymEngine::rcp_dynamic_cast<const SymEngine::Number>(
                    product->get_scalar()
                );
                auto factors = product->get_factors();
                auto chi = SymEngine::rcp_dynamic_cast<const TemporumOperator>(factors[0]);
                return chi->get_type() == TemporumType::Bra
                    ? std::make_tuple(
                          coef,
                          chi,
                          SymEngine::rcp_dynamic_cast<const TemporumOperator>(factors[1])
                      )
                    : std::make_tuple(
                          coef,
                          SymEngine::rcp_dynamic_cast<const TemporumOperator>(factors[1]),
                          chi
                      );
            }

            // Get frequency factor [sum(w_bra)-sum(w_ket)]/2 of a product
            inline SymEngine::RCP<const SymEngine::Number>
            get_frequency(const std::size_t index) const
            {
                auto term = get_braket_product(index);
                return SymEngine::divnum(
                    SymEngine::mulnum(
                        SymEngine::addnum(
                            std::get<1>(term)->get_frequency(),
                            std::get<2>(term)->get_frequency()
                        ),
                        std::get<0>(term)
                    ),
                    SymEngine::integer(-2)
                );
            }

            // Get derivatives on bra and ket of a product
            inline std::pair<SymEngine::multiset_basic, SymEngine::multiset_basic>
            get_derivatives(const std::size_t index) const
            {
                auto term = get_braket_product(index);
                return std::make_pair(
                    std::get<1>(term)->get_derivatives(),
                    std::get<2>(term)->get_derivatives()
                );
            }
    };

    // Helper function to make T matrix
    inline SymEngine::RCP<const TemporumOverlap> make_t_matrix(
        const PertDependency& dependencies = {}
    )
    {
        return SymEngine::make_rcp<const TemporumOverlap>(dependencies);
    }
}
