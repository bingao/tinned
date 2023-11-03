#include <cstddef>
#include <set>
#include <string>
#include <utility>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"

using namespace Tinned;

TEST_CASE("Test OneElecDensity and make_1el_density()", "[ElectronicState]")
{
    auto D_name = std::string("D");
    auto D = make_1el_density(D_name);
    REQUIRE(D->get_name() == D_name);
    auto derivative = D->get_derivative();
    REQUIRE(SymEngine::unified_eq(derivative, SymEngine::multiset_basic({})));

    auto el0 = make_perturbation(std::string("EL"), SymEngine::two);
    auto el1 = make_perturbation(
        std::string("EL"), SymEngine::two, std::set<std::size_t>({0, 1})
    );
    auto geo = make_perturbation(std::string("GEO"));
    auto Dp = SymEngine::rcp_dynamic_cast<const OneElecDensity>(
        (((D->diff(geo))->diff(el0))->diff(el1))->diff(el0)
    );
    REQUIRE(Dp->get_name() == D_name);
    derivative = Dp->get_derivative();
    REQUIRE(SymEngine::unified_eq(
        derivative, SymEngine::multiset_basic({el0, el0, el1, geo})
    ));
}

TEST_CASE("Test OneElecOperator and make_1el_operator()", "[OneElecOperator]")
{
    auto el = make_perturbation(std::string("EL"), SymEngine::two);
    auto geo = make_perturbation(std::string("GEO"));
    auto mag = make_perturbation(std::string("MAG"), SymEngine::one);
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto W_name = std::string("Omega");
    auto W = make_1el_operator(W_name, dependencies);
    REQUIRE(W->get_name() == W_name);
    REQUIRE(eq_dependency(dependencies, W->get_dependencies()));

    auto Wp = SymEngine::rcp_dynamic_cast<const OneElecOperator>(
        ((W->diff(geo))->diff(mag))->diff(geo)
    );
    REQUIRE(Wp->get_name() == W_name);
    auto derivative = Wp->get_derivative();
    REQUIRE(SymEngine::unified_eq(
        derivative, SymEngine::multiset_basic({geo, geo, mag})
    ));

    REQUIRE(SymEngine::eq(*Wp->diff(el), *make_zero_operator()));
}
