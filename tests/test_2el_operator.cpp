#include <map>
#include <string>
#include <utility>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_add.h>

#include "Tinned.hpp"

using namespace Tinned;

TEST_CASE("Test TwoElecOperator and make_2el_operator()", "[TwoElecOperator]")
{
    auto D_name = std::string("D");
    auto D = make_1el_density(D_name);
    auto el = make_perturbation(std::string("EL"), SymEngine::two);
    auto geo = make_perturbation(std::string("GEO"));
    auto mag = make_perturbation(std::string("MAG"), SymEngine::one);
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto G_name = std::string("G");
    auto G = make_2el_operator(G_name, D, dependencies);

    REQUIRE(G->get_name() == G_name);
    REQUIRE(SymEngine::eq(*D, *G->get_state()));
    REQUIRE(SymEngine::eq(*D, *G->get_args()[0]));
    REQUIRE(eq_dependency(dependencies, G->get_dependencies()));

    auto result = ((G->diff(el))->diff(geo))->diff(mag);
    // The above third order deritive should be
    // contr(g->diff(geo, mag), D->diff(el)) + contr(g->diff(geo), D->diff(el, mag))
    // + contr(g->diff(mag), D->diff(el, geo)) + contr(g, D->diff(el, geo, mag))
    REQUIRE(SymEngine::is_a_sub<const SymEngine::MatrixAdd>(*result));

    auto Df = SymEngine::rcp_dynamic_cast<const OneElecDensity>(D->diff(el));
    auto Dfg = SymEngine::rcp_dynamic_cast<const OneElecDensity>(Df->diff(geo));
    auto Dfb = SymEngine::rcp_dynamic_cast<const OneElecDensity>(Df->diff(mag));
    auto Dfgb = SymEngine::rcp_dynamic_cast<const OneElecDensity>(Dfg->diff(mag));
    auto found = std::map<SymEngine::RCP<const OneElecDensity>, bool, SymEngine::RCPBasicKeyLess>
    ({
        {Df, false}, {Dfg, false}, {Dfb, false}, {Dfgb, false}
    });
    auto args = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixAdd>(result)->get_args();
    for (auto& arg: args) {
        REQUIRE(SymEngine::is_a_sub<const TwoElecOperator>(*arg));
        auto Gp = SymEngine::rcp_dynamic_cast<const TwoElecOperator>(arg);
        REQUIRE(Gp->get_name() == G_name);
        auto Dp = Gp->get_state();
        REQUIRE(Dp->get_name() == D_name);
        auto derivatives = Gp->get_derivatives();
        if (SymEngine::eq(*Dp, *Df)) {
            REQUIRE(SymEngine::unified_eq(
                derivatives, SymEngine::multiset_basic({geo, mag})
            ));
            found[Df] = true;
        }
        else if (SymEngine::eq(*Dp, *Dfg)) {
            REQUIRE(SymEngine::unified_eq(
                derivatives, SymEngine::multiset_basic({mag})
            ));
            found[Dfg] = true;
        }
        else if (SymEngine::eq(*Dp, *Dfb)) {
            REQUIRE(SymEngine::unified_eq(
                derivatives, SymEngine::multiset_basic({geo})
            ));
            found[Dfb] = true;
        }
        else if (SymEngine::eq(*Dp, *Dfgb)) {
            REQUIRE(SymEngine::unified_eq(
                derivatives, SymEngine::multiset_basic({})
            ));
            found[Dfgb] = true;
        }
    }
    for (auto it = found.begin(); it != found.end(); ++it) REQUIRE(it->second);
}

TEST_CASE("Test TwoElecEnergy and make_2el_energy()", "[TwoElecEnergy]")
{
}
