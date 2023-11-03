#include <set>
#include <string>
#include <utility>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/integer.h>
#include <symengine/real_double.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"

using namespace Tinned;

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

TEST_CASE("Test ExchCorrEnergy and make_xc_energy()", "[ExchCorrEnergy]")
{
    auto D = make_1el_density(std::string("D"));
    auto a = make_perturbation(std::string("a"));
    auto b = make_perturbation(std::string("b"));
    auto c = make_perturbation(std::string("c"));
    auto d = make_perturbation(std::string("d"));
    auto Omega = make_1el_operator(
        std::string("Omega"),
        PertDependency({
            std::make_pair(a, 99),
            std::make_pair(b, 99),
            std::make_pair(c, 99),
            std::make_pair(d, 99)
        })
    );
    // We first do not consider derivatives of grid weights
    auto weight = make_nonel_function(std::string("weight"));
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
    // For unperturbed XC energy functional, there is only one unperturbed grid
    // weight, and order of the XC energy functional derivative is 0, and no
    // generalized density vectors to be contracted with
    auto terms = Exc->get_energy_terms();
    REQUIRE(terms.size() == 1);
    auto iter = terms.find(weight);
    REQUIRE(iter != terms.end());
    REQUIRE(iter->second.size() == 1);
    REQUIRE(iter->second.begin()->first == 0);
    REQUIRE(iter->second.begin()->second.is_null());

    // Tests from J. Chem. Phys. 140, 034103 (2014), equations (43)-(50).
    auto D_a = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(a));
    auto D_b = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(b));
    auto D_c = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(c));
    auto D_d = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(d));
    auto D_ab = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_a->diff(b));
    auto D_ac = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_a->diff(c));
    auto D_ad = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_a->diff(d));
    auto D_bc = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_b->diff(c));
    auto D_bd = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_b->diff(d));
    auto D_cd = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_c->diff(d));
    auto D_abc = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_ab->diff(c));
    auto D_abd = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_ab->diff(d));
    auto D_acd = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_ac->diff(d));
    auto D_bcd = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_bc->diff(d));
    auto D_abcd = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_abc->diff(d));
    auto Omega_a = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(a));
    auto Omega_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(b));
    auto Omega_c = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(c));
    auto Omega_d = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(d));
    auto Omega_ab = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_a->diff(b));
    auto Omega_ac = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_a->diff(c));
    auto Omega_ad = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_a->diff(d));
    auto Omega_bc = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_b->diff(c));
    auto Omega_bd = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_b->diff(d));
    auto Omega_cd = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_c->diff(d));
    auto Omega_abc = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_ab->diff(c));
    auto Omega_abd = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_ab->diff(d));
    auto Omega_acd = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_ac->diff(d));
    auto Omega_bcd = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_bc->diff(d));
    auto Omega_abcd = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_abc->diff(d));
    // The first order XC energy density derivative
    auto Exc_a = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc->diff(a));
    weights = Exc_a->get_weights();
    REQUIRE(weights.size() == 1);
    REQUIRE(SymEngine::eq(*(*weights.begin()), *weight));
    states = Exc_a->get_states();
    REQUIRE(states.size() == 2);
    REQUIRE(SymEngine::unified_eq(
        states,
        std::set<SymEngine::RCP<const ElectronicState>, SymEngine::RCPBasicKeyLess>({
          D, D_a
        })
    ));
    Omegas = Exc_a->get_overlap_distributions();
    REQUIRE(Omegas.size() == 2);
    REQUIRE(SymEngine::unified_eq(
        Omegas,
        std::set<SymEngine::RCP<const OneElecOperator>, SymEngine::RCPBasicKeyLess>({
          Omega, Omega_a
        })
    ));
    orders = Exc_a->get_exc_orders();
    REQUIRE(orders.size() == 1);
    REQUIRE(*orders.begin() == 1);
    terms = Exc_a->get_energy_terms();
    REQUIRE(terms.size() == 1);
    iter = terms.find(weight);
    REQUIRE(iter != terms.end());
    REQUIRE(iter->second.size() == 1);
    REQUIRE(iter->second.begin()->first == 1);
    REQUIRE(SymEngine::eq(
        *(iter->second.begin()->second),
        *SymEngine::add(
            SymEngine::mul(Omega_a, D),
            SymEngine::mul(Omega, D_a)
        )
    ));
}

TEST_CASE("Test ExchCorrPotential and make_xc_potential()", "[ExchCorrPotential]")
{
//    auto Fxc = SymEngine::make_rcp<const ExchCorrPotential>(
//        std::string("Fxc"), D
//    );
}
