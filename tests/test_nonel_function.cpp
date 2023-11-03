#include <string>
#include <utility>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"

using namespace Tinned;

TEST_CASE("Test NonElecFunction and make_nonel_function()", "[NonElecFunction]")
{
    auto el = make_perturbation(std::string("EL"), SymEngine::two);
    auto geo = make_perturbation(std::string("GEO"));
    auto mag = make_perturbation(std::string("MAG"), SymEngine::one);
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto h_name = std::string("hnuc");
    auto h = make_nonel_function(h_name, dependencies);
    REQUIRE(h->get_name() == h_name);
    REQUIRE(eq_dependency(dependencies, h->get_dependencies()));

    auto hp = SymEngine::rcp_dynamic_cast<const NonElecFunction>(
        ((h->diff(geo))->diff(mag))->diff(geo)
    );
    REQUIRE(hp->get_name() == h_name);
    auto derivative = hp->get_derivative();
    REQUIRE(SymEngine::unified_eq(
        derivative, SymEngine::multiset_basic({geo, geo, mag})
    ));

    REQUIRE(SymEngine::eq(*hp->diff(el), *SymEngine::zero));
}
