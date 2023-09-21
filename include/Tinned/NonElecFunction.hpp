/* tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of non-electron like functions.

   2023-09-08, Bin Gao:
   * first version
*/

#pragma once

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/functions.h>
#include <symengine/number.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

#include <symengine/eval.h>
#include <symengine/real_double.h>

// For example, the internuclear repulsion and nucleus interaction with
// external fields
class NonElecFunction: public SymEngine::FunctionWrapper
{
    public:
        OneElecOperator(const SymEngine::vec_basic &arg)
            : SymEngine::FunctionWrapper("OneElecOperator", arg) {}
        SymEngine::RCP<const SymEngine::Basic> create(const SymEngine::vec_basic &v) const
        {
            return SymEngine::make_rcp<OneElecOperator>(v);
        }
        SymEngine::RCP<const SymEngine::Number> eval(long bits) const
        {
            //FIXME: eval() will be performed for getting integrals or
            //expectation values with the knowledge of basis functions
            std::cout << "OneElecOperator::eval called\n";
    //RCP<const Symbol> x = symbol("x");
    //RCP<const Basic> e = add(one, make_rcp<MySin>(x));
    //RCP<const Basic> f;
    //f = e->subs({{x, integer(1)}});
    //double d = eval_double(*f);
            return SymEngine::real_double(3.14);
            //return real_double(::sin(eval_double(*get_vec()[0])));
        }
        SymEngine::RCP<const SymEngine::Basic> diff_impl(
            const SymEngine::RCP<const SymEngine::Symbol> &s
        ) const
        {
            std::cout << "OneElecOperator::diff_impl called\n";
            //FIXME: check if the operator depends on s by calling get_vec()
            return SymEngine::Derivative::create(rcp_from_this(), {s});
        }
};
