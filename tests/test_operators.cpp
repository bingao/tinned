#define CATCH_CONFIG_MAIN

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/number.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/real_double.h>
#include <symengine/complex.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_add.h>

#include "Tinned.hpp"

using namespace Tinned;

TEST_CASE("Test Perturbation and PertDependency", "[Perturbation]")
{
    auto el_name = std::string("EL");
    auto int_freq = SymEngine::zero;
    auto rat_freq = SymEngine::rational(1, 3);
    auto real_freq = SymEngine::real_double(0.5);
    auto cmplx_freq = SymEngine::Complex::from_two_nums(*rat_freq, *rat_freq);

    auto el0 = SymEngine::make_rcp<const Perturbation>(el_name, int_freq);
    auto el1 = SymEngine::make_rcp<const Perturbation>(el_name, rat_freq);
    auto el2 = SymEngine::make_rcp<const Perturbation>(el_name, real_freq);
    auto el3 = SymEngine::make_rcp<const Perturbation>(el_name, cmplx_freq);
    REQUIRE(el0->get_name() == el_name);
    REQUIRE(SymEngine::neq(*el0, *el1));
    REQUIRE(SymEngine::neq(*el0, *el2));
    REQUIRE(SymEngine::neq(*el0, *el3));
    REQUIRE(SymEngine::neq(*el1, *el2));
    REQUIRE(SymEngine::neq(*el1, *el3));
    REQUIRE(SymEngine::neq(*el2, *el3));

    auto components = std::set<std::size_t>({0, 1});
    auto el4 = SymEngine::make_rcp<const Perturbation>(el_name, int_freq, components);
    REQUIRE(SymEngine::neq(*el0, *el4));
    REQUIRE(SymEngine::neq(*el1, *el4));
    REQUIRE(SymEngine::neq(*el2, *el4));

    auto geo_name = std::string("GEO");
    auto geo0 = SymEngine::make_rcp<const Perturbation>(geo_name, int_freq);
    REQUIRE(SymEngine::neq(*el0, *geo0));
    REQUIRE(SymEngine::neq(*el1, *geo0));
    REQUIRE(SymEngine::neq(*el2, *geo0));
    REQUIRE(SymEngine::neq(*el4, *geo0));

    auto geo1 = SymEngine::make_rcp<const Perturbation>(geo_name, int_freq);
    REQUIRE(SymEngine::eq(*geo0, *geo1));

    auto freq = el0->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *int_freq));
    freq = el1->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *rat_freq));
    freq = el2->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *real_freq));
    freq = el3->get_frequency();
    REQUIRE(SymEngine::eq(*freq, *cmplx_freq));

    REQUIRE(el0->get_components() == std::set<std::size_t>({}));
    REQUIRE(el4->get_components() == components);

    auto el_pert = PertDependency({
        std::make_pair(el0, 0),
        std::make_pair(el1, 1),
        std::make_pair(el2, 2),
        std::make_pair(el3, 3),
        std::make_pair(el4, 4)
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

    auto s = SymEngine::make_rcp<const SymEngine::Symbol>("s");
    REQUIRE(find_dependency(el_geo_pert, s) == 0);
    REQUIRE(find_dependency(el_geo_pert, el0) == 0);
    REQUIRE(find_dependency(el_geo_pert, el1) == 1);
    REQUIRE(find_dependency(el_geo_pert, el3) == 3);
    REQUIRE(find_dependency(el_geo_pert, geo0) == 99);
}

TEST_CASE("Test OneElecDensity", "[ElectronicState]")
{
    auto D_name = std::string("D");
    auto D = SymEngine::make_rcp<const OneElecDensity>(D_name);
    REQUIRE(D->get_name() == D_name);
    auto derivative = D->get_derivative();
    REQUIRE(SymEngine::unified_eq(derivative, SymEngine::multiset_basic({})));

    auto el0 = SymEngine::make_rcp<const Perturbation>(
        std::string("EL"), SymEngine::two
    );
    auto el1 = SymEngine::make_rcp<const Perturbation>(
        std::string("EL"), SymEngine::two, std::set<std::size_t>({0, 1})
    );
    auto geo = SymEngine::make_rcp<const Perturbation>(
        std::string("GEO"), SymEngine::zero
    );
    auto Dp = SymEngine::rcp_dynamic_cast<const OneElecDensity>(
        (((D->diff(geo))->diff(el0))->diff(el1))->diff(el0)
    );
    REQUIRE(Dp->get_name() == D_name);
    derivative = Dp->get_derivative();
    REQUIRE(SymEngine::unified_eq(
        derivative, SymEngine::multiset_basic({el0, el0, el1, geo})
    ));
}

TEST_CASE("Test OneElecOperator and make_overlap_distribution", "[OneElecOperator]")
{
    auto el = SymEngine::make_rcp<const Perturbation>(
        std::string("EL"), SymEngine::two
    );
    auto geo = SymEngine::make_rcp<const Perturbation>(
        std::string("GEO"), SymEngine::zero
    );
    auto mag = SymEngine::make_rcp<const Perturbation>(
        std::string("MAG"), SymEngine::one
    );
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto W_name = std::string("Omega");
    auto W = make_overlap_distribution(W_name, dependencies);
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

TEST_CASE("Test TwoElecOperator", "[TwoElecOperator]")
{
    auto D_name = std::string("D");
    auto D = SymEngine::make_rcp<const OneElecDensity>(D_name);
    auto el = SymEngine::make_rcp<const Perturbation>(
        std::string("EL"), SymEngine::two
    );
    auto geo = SymEngine::make_rcp<const Perturbation>(
        std::string("GEO"), SymEngine::zero
    );
    auto mag = SymEngine::make_rcp<const Perturbation>(
        std::string("MAG"), SymEngine::one
    );
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto G_name = std::string("G");
    auto G = SymEngine::make_rcp<const TwoElecOperator>(G_name, D, dependencies);

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
    auto args
        = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixAdd>(result)->get_args();
    for (auto& arg: args) {
        REQUIRE(SymEngine::is_a_sub<const TwoElecOperator>(*arg));
        auto Gp = SymEngine::rcp_dynamic_cast<const TwoElecOperator>(arg);
        REQUIRE(Gp->get_name() == G_name);
        auto Dp = Gp->get_state();
        REQUIRE(Dp->get_name() == D_name);
        auto derivative = Gp->get_derivative();
        if (SymEngine::eq(*Dp, *Df)) {
            REQUIRE(SymEngine::unified_eq(
                derivative, SymEngine::multiset_basic({geo, mag})
            ));
            found[Df] = true;
        }
        else if (SymEngine::eq(*Dp, *Dfg)) {
            REQUIRE(SymEngine::unified_eq(
                derivative, SymEngine::multiset_basic({mag})
            ));
            found[Dfg] = true;
        }
        else if (SymEngine::eq(*Dp, *Dfb)) {
            REQUIRE(SymEngine::unified_eq(
                derivative, SymEngine::multiset_basic({geo})
            ));
            found[Dfb] = true;
        }
        else if (SymEngine::eq(*Dp, *Dfgb)) {
            REQUIRE(SymEngine::unified_eq(
                derivative, SymEngine::multiset_basic({})
            ));
            found[Dfgb] = true;
        }
    }
    for (auto it = found.begin(); it != found.end(); ++it) REQUIRE(it->second);
}

TEST_CASE("Test CompositeFunction and ExchCorrEnergy", "[ExchCorrEnergy]")
{
//    auto Exc = SymEngine::make_rcp<const ExchCorrEnergy>(
//        std::string("Exc"), D
//    );
}

TEST_CASE("Test ExchCorrPotential", "[ExchCorrPotential]")
{
//    auto Fxc = SymEngine::make_rcp<const ExchCorrPotential>(
//        std::string("Fxc"), D
//    );
}

TEST_CASE("Test NonElecFunction", "[NonElecFunction]")
{
    auto el = SymEngine::make_rcp<const Perturbation>(
        std::string("EL"), SymEngine::two
    );
    auto geo = SymEngine::make_rcp<const Perturbation>(
        std::string("GEO"), SymEngine::zero
    );
    auto mag = SymEngine::make_rcp<const Perturbation>(
        std::string("MAG"), SymEngine::one
    );
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto h_name = std::string("hnuc");
    auto h = SymEngine::make_rcp<const NonElecFunction>(h_name, dependencies);
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

TEST_CASE("Test TemporumOperator", "[TemporumOperator]")
{
    auto D = SymEngine::make_rcp<const OneElecDensity>(std::string("D"));
    auto Dbra = SymEngine::make_rcp<const TemporumOperator>(D, TemporumType::Bra);
    auto Dket = SymEngine::make_rcp<const TemporumOperator>(D, TemporumType::Ket);
    REQUIRE(SymEngine::eq(*D, *Dbra->get_args()[0]));
    REQUIRE(SymEngine::eq(*D, *Dket->get_args()[0]));
    REQUIRE(SymEngine::eq(*D, *Dbra->get_target()));
    REQUIRE(SymEngine::eq(*D, *Dket->get_target()));
    REQUIRE(Dbra->get_type() == TemporumType::Bra);
    REQUIRE(Dket->get_type() == TemporumType::Ket);
    REQUIRE(SymEngine::eq(*Dbra->get_frequency(), *SymEngine::zero));
    REQUIRE(SymEngine::eq(*Dket->get_frequency(), *SymEngine::zero));

    auto el = SymEngine::make_rcp<const Perturbation>(
        std::string("EL"), SymEngine::real_double(0.5)
    );
    auto geo = SymEngine::make_rcp<const Perturbation>(
        std::string("GEO"), SymEngine::real_double(0.0)
    );
    auto mag = SymEngine::make_rcp<const Perturbation>(
        std::string("MAG"), SymEngine::real_double(1.5)
    );
    auto Dp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        ((Dbra->diff(el))->diff(geo))->diff(mag)
    );
    auto sum_freq = SymEngine::addnum(
        SymEngine::addnum(el->get_frequency(), geo->get_frequency()),
        mag->get_frequency()
    );
    REQUIRE(SymEngine::eq(
        *Dp->get_frequency(),
        *SymEngine::subnum(SymEngine::real_double(0), sum_freq)
    ));
    REQUIRE(SymEngine::unified_eq(
        Dp->get_derivative(), SymEngine::multiset_basic({el, geo, mag})
    ));
    Dp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        ((Dket->diff(el))->diff(geo))->diff(mag)
    );
    REQUIRE(SymEngine::eq(*Dp->get_frequency(), *sum_freq));
    REQUIRE(SymEngine::unified_eq(
        Dp->get_derivative(), SymEngine::multiset_basic({el, geo, mag})
    ));

    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto S = SymEngine::make_rcp<const OneElecOperator>(
        std::string("S"), dependencies
    );
    auto St = SymEngine::make_rcp<const TemporumOperator>(S);
    auto Sp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        (St->diff(geo))->diff(mag)
    );
    REQUIRE(SymEngine::unified_eq(
        Sp->get_derivative(), SymEngine::multiset_basic({geo, mag})
    ));
    REQUIRE(SymEngine::eq(*Sp->diff(el), *make_zero_operator()));

    auto h = SymEngine::make_rcp<const NonElecFunction>(
        std::string("hnuc"), dependencies
    );
    auto ht = SymEngine::make_rcp<const TemporumOperator>(h);
    auto hp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        (ht->diff(geo))->diff(mag)
    );
    REQUIRE(SymEngine::unified_eq(
        hp->get_derivative(), SymEngine::multiset_basic({geo, mag})
    ));
    REQUIRE(SymEngine::eq(*hp->diff(el), *SymEngine::zero));
}

TEST_CASE("Test TemporumOverlap", "[TemporumOverlap]")
{
    auto el_freq = SymEngine::real_double(0.5);
    auto el = SymEngine::make_rcp<const Perturbation>(std::string("EL"), el_freq);
    auto geo_freq = SymEngine::real_double(1.5);
    auto geo = SymEngine::make_rcp<const Perturbation>(std::string("GEO"), geo_freq);
    auto mag_freq = SymEngine::real_double(2.5);
    auto mag = SymEngine::make_rcp<const Perturbation>(std::string("MAG"), mag_freq);
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto T = SymEngine::make_rcp<const TemporumOverlap>(dependencies);
    REQUIRE(eq_dependency(dependencies, T->get_dependencies()));

    // = wg*S^{gg|0} + 0*2*S^{g|g} - wg*S^{0|gg}
    auto Tgg = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(
        (T->diff(geo))->diff(geo)
    );
    auto found = std::map<std::string, bool>({
        {std::string("gg|0"), false},
        {std::string("g|g"), false},
        {std::string("0|gg"), false}
    });
    REQUIRE(Tgg->size() == 3);
    for (std::size_t i = 0; i < Tgg->size(); ++i) {
        auto derivative = Tgg->get_derivative(i);
        if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo, geo})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({})
            )
        ) {
            REQUIRE(SymEngine::eq(*Tgg->get_frequency(i), *geo_freq));
            found[std::string("gg|0")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo})
            )
        ) {
            REQUIRE(Tgg->get_frequency(i)->is_zero());
            found[std::string("g|g")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo, geo})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tgg->get_frequency(i),
                *SymEngine::subnum(SymEngine::real_double(0.0), geo_freq)
            ));
            found[std::string("0|gg")] = true;
        }
    }
    for (auto it = found.begin(); it != found.end(); ++it) REQUIRE(it->second);

    // = (wg+wb/2)*S^{ggb|0} + (wg-wb/2)*S^{gg|b} + wb*S^{gb|g}
    // - wb*S^{g|gb} + (wb/2-wg)*S^{b|gg} - (wg+wb/2)*S^{0|ggb}
    auto Tggb = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(Tgg->diff(mag));
    found = std::map<std::string, bool>({
        {std::string("ggb|0"), false},
        {std::string("gg|b"), false},
        {std::string("gb|g"), false},
        {std::string("g|gb"), false},
        {std::string("b|gg"), false},
        {std::string("0|ggb"), false}
    });
    REQUIRE(Tggb->size() == 6);
    for (std::size_t i = 0; i < Tggb->size(); ++i) {
        auto derivative = Tggb->get_derivative(i);
        if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo, geo, mag})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::addnum(
                    geo_freq,
                    SymEngine::mulnum(mag_freq, SymEngine::real_double(0.5))
                )
            ));
            found[std::string("ggb|0")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo, geo})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({mag})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::subnum(
                    geo_freq,
                    SymEngine::mulnum(mag_freq, SymEngine::real_double(0.5))
                )
            ));
            found[std::string("gg|b")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo, mag})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo})
            )
        ) {
            REQUIRE(SymEngine::eq(*Tggb->get_frequency(i), *mag_freq));
            found[std::string("gb|g")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo, mag})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::subnum(SymEngine::real_double(0.0), mag_freq)
            ));
            found[std::string("g|gb")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({mag})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo, geo})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::subnum(
                    SymEngine::mulnum(mag_freq, SymEngine::real_double(0.5)),
                    geo_freq
                )
            ));
            found[std::string("b|gg")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo, geo, mag})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::subnum(
                    SymEngine::mulnum(mag_freq, SymEngine::real_double(-0.5)),
                    geo_freq
                )
            ));
            found[std::string("0|ggb")] = true;
        }
    }
    for (auto it = found.begin(); it != found.end(); ++it) REQUIRE(it->second);

    REQUIRE(SymEngine::eq(*T->diff(el), *make_zero_operator()));
    REQUIRE(SymEngine::eq(*Tgg->diff(el), *make_zero_operator()));
    REQUIRE(SymEngine::eq(*Tggb->diff(el), *make_zero_operator()));
}
