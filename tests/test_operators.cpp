#define CATCH_CONFIG_MAIN

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/number.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/real_double.h>
#include <symengine/complex.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_add.h>

#include "Tinned.hpp"

using namespace Tinned;

TEST_CASE("Test Perturbation, make_perturbation() and PertDependency", "[Perturbation]")
{
    auto el_name = std::string("EL");
    auto int_freq = SymEngine::zero;
    auto rat_freq = SymEngine::rational(1, 3);
    auto real_freq = SymEngine::real_double(0.5);
    auto cmplx_freq = SymEngine::Complex::from_two_nums(*rat_freq, *rat_freq);

    auto el0 = make_perturbation(el_name, int_freq);
    auto el1 = make_perturbation(el_name, rat_freq);
    auto el2 = make_perturbation(el_name, real_freq);
    auto el3 = make_perturbation(el_name, cmplx_freq);
    REQUIRE(el0->get_name() == el_name);
    REQUIRE(SymEngine::neq(*el0, *el1));
    REQUIRE(SymEngine::neq(*el0, *el2));
    REQUIRE(SymEngine::neq(*el0, *el3));
    REQUIRE(SymEngine::neq(*el1, *el2));
    REQUIRE(SymEngine::neq(*el1, *el3));
    REQUIRE(SymEngine::neq(*el2, *el3));

    auto components = std::set<std::size_t>({0, 1});
    auto el4 = make_perturbation(el_name, int_freq, components);
    REQUIRE(SymEngine::neq(*el0, *el4));
    REQUIRE(SymEngine::neq(*el1, *el4));
    REQUIRE(SymEngine::neq(*el2, *el4));

    auto geo_name = std::string("GEO");
    auto geo0 = make_perturbation(geo_name, int_freq);
    REQUIRE(SymEngine::neq(*el0, *geo0));
    REQUIRE(SymEngine::neq(*el1, *geo0));
    REQUIRE(SymEngine::neq(*el2, *geo0));
    REQUIRE(SymEngine::neq(*el4, *geo0));

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

    auto s = SymEngine::symbol("s");
    REQUIRE(find_dependency(el_geo_pert, s) == 0);
    REQUIRE(find_dependency(el_geo_pert, el0) == 0);
    REQUIRE(find_dependency(el_geo_pert, el1) == 1);
    REQUIRE(find_dependency(el_geo_pert, el3) == 3);
    REQUIRE(find_dependency(el_geo_pert, geo0) == 99);
}

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

TEST_CASE("Test CompositeFunction", "[CompositeFunction]")
{
    auto el_freq = SymEngine::real_double(0.5);
    auto el = make_perturbation(std::string("EL"), el_freq);
    auto geo_freq = SymEngine::real_double(1.5);
    auto geo = make_perturbation(std::string("GEO"), geo_freq);
    auto mag_freq = SymEngine::real_double(2.5);
    auto mag = make_perturbation(std::string("MAG"), mag_freq);
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto inner_name = std::string("inner");
    auto inner = make_nonel_function(inner_name, dependencies);
    auto fun_name = std::string("outer");
    auto fun = SymEngine::make_rcp<const CompositeFunction>(fun_name, inner);
    REQUIRE(fun->get_name() == fun_name);
    REQUIRE(SymEngine::eq(*fun->get_inner(), *inner));

    // = outer^{(3)}*inner^{g}*inner^{g}*inner^{b} + 2*outer^{(2)}*inner^{gb}*inner^{g}
    // + outer^{(2)}*inner^{gg}*inner^{b} + outer^{(1)}*inner^{ggb}
    auto fun_ggm = ((fun->diff(geo))->diff(geo))->diff(mag);
    auto fun_1 = SymEngine::make_rcp<const CompositeFunction>(fun_name, inner, 1);
    auto fun_2 = SymEngine::make_rcp<const CompositeFunction>(fun_name, inner, 2);
    auto fun_3 = SymEngine::make_rcp<const CompositeFunction>(fun_name, inner, 3);
    auto inner_g = inner->diff(geo);
    auto inner_b = inner->diff(mag);
    auto inner_gg = inner_g->diff(geo);
    auto inner_gb = inner_g->diff(mag);
    auto inner_ggb = inner_gg->diff(mag);
    REQUIRE(SymEngine::eq(
        *fun_ggm,
        *SymEngine::add(SymEngine::vec_basic({
            SymEngine::mul(SymEngine::vec_basic({
                fun_3, inner_g, inner_g, inner_b
            })),
            SymEngine::mul(SymEngine::vec_basic({
                SymEngine::integer(2), fun_2, inner_gb, inner_g
            })),
            SymEngine::mul(SymEngine::vec_basic({
                fun_2, inner_gg, inner_b
            })),
            SymEngine::mul(SymEngine::vec_basic({
                fun_1, inner_ggb
            }))
        }))
    ));

    // Test from https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula#Example
    auto fun_gggg = (((fun->diff(geo))->diff(geo))->diff(geo))->diff(geo);
    auto fun_4 = SymEngine::make_rcp<const CompositeFunction>(fun_name, inner, 4);
    auto inner_ggg = inner_gg->diff(geo);
    auto inner_gggg = inner_ggg->diff(geo);
    REQUIRE(SymEngine::eq(
        *fun_gggg,
        *SymEngine::add(SymEngine::vec_basic({
            SymEngine::mul(SymEngine::vec_basic({
                fun_4, inner_g, inner_g, inner_g, inner_g
            })),
            SymEngine::mul(SymEngine::vec_basic({
                SymEngine::integer(6), fun_3, inner_gg, inner_g, inner_g
            })),
            SymEngine::mul(SymEngine::vec_basic({
                SymEngine::integer(3), fun_2, inner_gg, inner_gg
            })),
            SymEngine::mul(SymEngine::vec_basic({
                SymEngine::integer(4), fun_2, inner_ggg, inner_g
            })),
            SymEngine::mul(SymEngine::vec_basic({fun_1, inner_gggg}))
        }))
    ));

    REQUIRE(SymEngine::eq(*fun->diff(el), *SymEngine::zero));
    REQUIRE(SymEngine::eq(*fun_ggm->diff(el), *SymEngine::zero));
    REQUIRE(SymEngine::eq(*fun_gggg->diff(el), *SymEngine::zero));
}

#include <iostream>

TEST_CASE("Test ExchCorrEnergy and make_xc_energy()", "[ExchCorrEnergy]")
{
    auto D = make_1el_density(std::string("D"));
    auto el = make_perturbation(std::string("EL"), SymEngine::two);
    auto geo = make_perturbation(std::string("GEO"));
    auto mag = make_perturbation(std::string("MAG"), SymEngine::one);
    auto Omega = make_1el_operator(
        std::string("Omega"),
        PertDependency({
            std::make_pair(geo, 99), std::make_pair(mag, 99)
        })
    );
    auto weight = make_nonel_function(
        std::string("weight"),
        PertDependency({std::make_pair(geo, 99)})
    );
    auto Exc_name = std::string("GGA");
    auto Exc = make_xc_energy(Exc_name, D, Omega, weight);

    REQUIRE(Exc->get_name() == Exc_name);
    REQUIRE(SymEngine::eq(*Exc->get_weight(), *weight));
    REQUIRE(SymEngine::eq(*Exc->get_state(), *D));
    REQUIRE(SymEngine::eq(*Exc->get_overlap_distribution(), *Omega));
    auto weights = Exc->get_weights();
    REQUIRE(weights.size() == 1);
    REQUIRE(SymEngine::eq(*(*weights.begin()), *weight));
    auto states = Exc->get_states();
    REQUIRE(states.size() == 1);
    REQUIRE(SymEngine::eq(*(*states.begin()), *D));
    auto Omegas = Exc->get_overlap_distributions();
    REQUIRE(Omegas.size() == 1);
    REQUIRE(SymEngine::eq(*(*Omegas.begin()), *Omega));
    auto orders = Exc->get_exc_orders();
    REQUIRE(orders.size() == 1);
    REQUIRE(*orders.begin() == 0);

    // Tests from J. Chem. Phys. 140, 034103 (2014)
    auto Exc_f = Exc->diff(el);
    auto Exc_fb = Exc_f->diff(mag);
    auto Exc_fbb = Exc_fb->diff(mag);
std::cout << "Exc: " << stringify(Exc) << "\n\n";
std::cout << "Exc->diff(el): " << stringify(Exc_f) << "\n\n";
std::cout << "(Exc->diff(el))->diff(mag): " << stringify(Exc_fb) << "\n\n";
std::cout << "((Exc->diff(el))->diff(mag))->diff(mag): " << stringify(Exc_fbb) << "\n\n";
    //get_energy()
    //get_weights()
    //get_states()
    //get_overlap_distributions()
    //get_exc_orders()
}

TEST_CASE("Test ExchCorrPotential and make_xc_potential()", "[ExchCorrPotential]")
{
//    auto Fxc = SymEngine::make_rcp<const ExchCorrPotential>(
//        std::string("Fxc"), D
//    );
}

TEST_CASE("Test TemporumOperator and make_dt_operator()", "[TemporumOperator]")
{
    auto D = make_1el_density(std::string("D"));
    auto Dbra = make_dt_operator(D, TemporumType::Bra);
    auto Dket = make_dt_operator(D, TemporumType::Ket);
    REQUIRE(SymEngine::eq(*D, *Dbra->get_args()[0]));
    REQUIRE(SymEngine::eq(*D, *Dket->get_args()[0]));
    REQUIRE(SymEngine::eq(*D, *Dbra->get_target()));
    REQUIRE(SymEngine::eq(*D, *Dket->get_target()));
    REQUIRE(Dbra->get_type() == TemporumType::Bra);
    REQUIRE(Dket->get_type() == TemporumType::Ket);
    REQUIRE(SymEngine::eq(*Dbra->get_frequency(), *SymEngine::zero));
    REQUIRE(SymEngine::eq(*Dket->get_frequency(), *SymEngine::zero));

    auto el = make_perturbation(std::string("EL"), SymEngine::real_double(0.5));
    auto geo = make_perturbation(std::string("GEO"), SymEngine::real_double(0.0));
    auto mag = make_perturbation(std::string("MAG"), SymEngine::real_double(1.5));
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
    auto S = make_1el_operator(std::string("S"), dependencies);
    auto St = make_dt_operator(S);
    auto Sp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        (St->diff(geo))->diff(mag)
    );
    REQUIRE(SymEngine::unified_eq(
        Sp->get_derivative(), SymEngine::multiset_basic({geo, mag})
    ));
    REQUIRE(SymEngine::eq(*Sp->diff(el), *make_zero_operator()));

    auto h = make_nonel_function(std::string("hnuc"), dependencies);
    auto ht = make_dt_operator(h);
    auto hp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        (ht->diff(geo))->diff(mag)
    );
    REQUIRE(SymEngine::unified_eq(
        hp->get_derivative(), SymEngine::multiset_basic({geo, mag})
    ));
    REQUIRE(SymEngine::eq(*hp->diff(el), *SymEngine::zero));
}

TEST_CASE("Test TemporumOverlap and make_t_matrix()", "[TemporumOverlap]")
{
    auto el_freq = SymEngine::real_double(0.5);
    auto el = make_perturbation(std::string("EL"), el_freq);
    auto geo_freq = SymEngine::real_double(1.5);
    auto geo = make_perturbation(std::string("GEO"), geo_freq);
    auto mag_freq = SymEngine::real_double(2.5);
    auto mag = make_perturbation(std::string("MAG"), mag_freq);
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto T = make_t_matrix(dependencies);
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
