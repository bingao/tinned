#define CATCH_CONFIG_MAIN

#include <cstddef>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/real_double.h>
#include <symengine/complex.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"

using namespace Tinned;

TEST_CASE("Test Perturbation, make_perturbation() and PertDependency", "[Perturbation]")
{
    auto el_name = std::string("EL");
    auto int_freq = SymEngine::zero;
    auto rat_freq = SymEngine::rational(1, 3);
    auto real_freq = SymEngine::real_double(0.5);
    auto cmplx_freq = SymEngine::Complex::from_two_nums(*rat_freq, *rat_freq);
    auto symbol_freq = SymEngine::symbol("omega");
    auto sum_freq = SymEngine::add(rat_freq, symbol_freq);

    auto el0 = make_perturbation(el_name, int_freq);
    auto el1 = make_perturbation(el_name, rat_freq);
    auto el2 = make_perturbation(el_name, real_freq);
    auto el3 = make_perturbation(el_name, cmplx_freq);
    auto el4 = make_perturbation(el_name, symbol_freq);
    auto el5 = make_perturbation(el_name, sum_freq);
    REQUIRE(el0->get_name() == el_name);
    REQUIRE(SymEngine::neq(*el0, *el1));
    REQUIRE(SymEngine::neq(*el0, *el2));
    REQUIRE(SymEngine::neq(*el0, *el3));
    REQUIRE(SymEngine::neq(*el0, *el4));
    REQUIRE(SymEngine::neq(*el0, *el5));
    REQUIRE(SymEngine::neq(*el1, *el2));
    REQUIRE(SymEngine::neq(*el1, *el3));
    REQUIRE(SymEngine::neq(*el1, *el4));
    REQUIRE(SymEngine::neq(*el1, *el5));
    REQUIRE(SymEngine::neq(*el2, *el3));
    REQUIRE(SymEngine::neq(*el2, *el4));
    REQUIRE(SymEngine::neq(*el2, *el5));
    REQUIRE(SymEngine::neq(*el3, *el4));
    REQUIRE(SymEngine::neq(*el3, *el5));
    REQUIRE(SymEngine::neq(*el4, *el5));

    auto components = std::set<std::size_t>({0, 1});
    auto el6 = make_perturbation(el_name, int_freq, components);
    REQUIRE(SymEngine::neq(*el0, *el6));
    REQUIRE(SymEngine::neq(*el1, *el6));
    REQUIRE(SymEngine::neq(*el2, *el6));
    REQUIRE(SymEngine::neq(*el3, *el6));
    REQUIRE(SymEngine::neq(*el4, *el6));
    REQUIRE(SymEngine::neq(*el5, *el6));

    auto geo_name = std::string("GEO");
    auto geo0 = make_perturbation(geo_name, int_freq);
    REQUIRE(SymEngine::neq(*el0, *geo0));
    REQUIRE(SymEngine::neq(*el1, *geo0));
    REQUIRE(SymEngine::neq(*el2, *geo0));
    REQUIRE(SymEngine::neq(*el3, *geo0));
    REQUIRE(SymEngine::neq(*el4, *geo0));
    REQUIRE(SymEngine::neq(*el5, *geo0));
    REQUIRE(SymEngine::neq(*el6, *geo0));

    auto geo1 = make_perturbation(geo_name, int_freq);
    REQUIRE(SymEngine::eq(*geo0, *geo1));

    auto freq = el0->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *int_freq));
    freq = el1->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *rat_freq));
    freq = el2->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *real_freq));
    freq = el3->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *cmplx_freq));
    freq = el4->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *symbol_freq));
    freq = el5->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *sum_freq));

    REQUIRE(el0->get_components() == std::set<std::size_t>({}));
    REQUIRE(el6->get_components() == components);

    auto el_pert = PertDependency({
        std::make_pair(el0, 0),
        std::make_pair(el1, 1),
        std::make_pair(el2, 2),
        std::make_pair(el3, 3),
        std::make_pair(el4, 4),
        std::make_pair(el5, 5),
        std::make_pair(el6, 6)
    });
    auto geo_pert = PertDependency({
        std::make_pair(geo0, 99)
    });
    auto el_geo_pert = PertDependency({
        std::make_pair(el0, 0),
        std::make_pair(el1, 1),
        std::make_pair(el3, 3),
        std::make_pair(geo0, 99)
    });

    REQUIRE(eq_dependency(el_pert, el_pert));
    REQUIRE(!eq_dependency(el_pert, geo_pert));
    REQUIRE(!eq_dependency(el_pert, el_geo_pert));

    auto s = SymEngine::symbol("s");
    REQUIRE(get_diff_order(s, el_geo_pert) == 0);
    REQUIRE(get_diff_order(el0, el_geo_pert) == 0);
    REQUIRE(get_diff_order(el1, el_geo_pert) == 1);
    REQUIRE(get_diff_order(el3, el_geo_pert) == 3);
    REQUIRE(get_diff_order(geo0, el_geo_pert) == 99);

    REQUIRE(!is_zero_derivative(SymEngine::multiset_basic({}), el_pert));
    REQUIRE(is_zero_derivative(SymEngine::multiset_basic({geo0}), el_pert));
    REQUIRE(is_zero_derivative(SymEngine::multiset_basic({el1, el1}), el_pert));
    REQUIRE(!is_zero_derivative(SymEngine::multiset_basic({el1, el3, el4, el6}), el_pert));
}
