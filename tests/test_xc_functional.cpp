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

    // Lambda expressions to make contractions of (un)perturbed generalized
    // density vectors and the XC functional derivative vectors
    //
    // For unperturbed XC energy functional, there is only one unperturbed grid
    // weight, and order of the XC energy functional derivative is 0, and no
    // generalized density vectors to be contracted with
    auto make_unperturbed_contraction = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const OneElecOperator>& Omega
    ) {
        return ExcDensityContractionMap({
            {make_exc_density(D, Omega, 0), SymEngine::RCP<const SymEngine::Basic>()}
        });
    };
    // Equation (47), J. Chem. Phys. 140, 034103 (2014)
    auto make_first_order_density = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a
    ) {
        return SymEngine::add(
            make_density_vector(D, Omega_a),
            make_density_vector(D_a, Omega)
        );
    };
    // Equation (43), J. Chem. Phys. 140, 034103 (2014)
    auto make_first_order_contraction = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a
    ) {
        return ExcDensityContractionMap({
            {
                make_exc_density(D, Omega, 1),
                make_first_order_density(D, D_a, Omega, Omega_a)
            }
        });
    };
    // Equation (48), J. Chem. Phys. 140, 034103 (2014)
    auto make_second_order_density = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const ElectronicState>& D_b,
        const SymEngine::RCP<const ElectronicState>& D_ab,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a,
        const SymEngine::RCP<const OneElecOperator>& Omega_b,
        const SymEngine::RCP<const OneElecOperator>& Omega_ab
    ) {
        return SymEngine::add({
            make_density_vector(D, Omega_ab),
            make_density_vector(D_a, Omega_b),
            make_density_vector(D_b, Omega_a),
            make_density_vector(D_ab, Omega)
        });
    };
    // Equation (44), J. Chem. Phys. 140, 034103 (2014)
    auto make_second_order_contraction = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const ElectronicState>& D_b,
        const SymEngine::RCP<const ElectronicState>& D_ab,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a,
        const SymEngine::RCP<const OneElecOperator>& Omega_b,
        const SymEngine::RCP<const OneElecOperator>& Omega_ab
    ) {
        return ExcDensityContractionMap({
            {
                make_exc_density(D, Omega, 1),
                make_second_order_density(
                    D, D_a, D_b, D_ab, Omega, Omega_a, Omega_b, Omega_ab
                )
            },
            {
                make_exc_density(D, Omega, 2),
                SymEngine::expand(SymEngine::mul(
                    make_first_order_density(D, D_a, Omega, Omega_a),
                    make_first_order_density(D, D_b, Omega, Omega_b)
                ))
            }
        });
    };
    // Equation (49), J. Chem. Phys. 140, 034103 (2014)
    auto make_third_order_density = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const ElectronicState>& D_b,
        const SymEngine::RCP<const ElectronicState>& D_c,
        const SymEngine::RCP<const ElectronicState>& D_ab,
        const SymEngine::RCP<const ElectronicState>& D_ac,
        const SymEngine::RCP<const ElectronicState>& D_bc,
        const SymEngine::RCP<const ElectronicState>& D_abc,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a,
        const SymEngine::RCP<const OneElecOperator>& Omega_b,
        const SymEngine::RCP<const OneElecOperator>& Omega_c,
        const SymEngine::RCP<const OneElecOperator>& Omega_ab,
        const SymEngine::RCP<const OneElecOperator>& Omega_ac,
        const SymEngine::RCP<const OneElecOperator>& Omega_bc,
        const SymEngine::RCP<const OneElecOperator>& Omega_abc
    ) {
        return SymEngine::add({
            make_density_vector(D, Omega_abc),
            make_density_vector(D_a, Omega_bc),
            make_density_vector(D_b, Omega_ac),
            make_density_vector(D_c, Omega_ab),
            make_density_vector(D_ab, Omega_c),
            make_density_vector(D_ac, Omega_b),
            make_density_vector(D_bc, Omega_a),
            make_density_vector(D_abc, Omega)
        });
    };
    // Equation (45), J. Chem. Phys. 140, 034103 (2014)
    auto make_third_order_contraction = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const ElectronicState>& D_b,
        const SymEngine::RCP<const ElectronicState>& D_c,
        const SymEngine::RCP<const ElectronicState>& D_ab,
        const SymEngine::RCP<const ElectronicState>& D_ac,
        const SymEngine::RCP<const ElectronicState>& D_bc,
        const SymEngine::RCP<const ElectronicState>& D_abc,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a,
        const SymEngine::RCP<const OneElecOperator>& Omega_b,
        const SymEngine::RCP<const OneElecOperator>& Omega_c,
        const SymEngine::RCP<const OneElecOperator>& Omega_ab,
        const SymEngine::RCP<const OneElecOperator>& Omega_ac,
        const SymEngine::RCP<const OneElecOperator>& Omega_bc,
        const SymEngine::RCP<const OneElecOperator>& Omega_abc
    ) {
        return ExcDensityContractionMap({
            {
                make_exc_density(D, Omega, 1),
                make_third_order_density(
                    D, D_a, D_b, D_c, D_ab, D_ac, D_bc, D_abc,
                    Omega,
                    Omega_a, Omega_b, Omega_c,
                    Omega_ab, Omega_ac, Omega_bc,
                    Omega_abc
                )
            },
            {
                make_exc_density(D, Omega, 2),
                SymEngine::expand(SymEngine::add({
                    SymEngine::mul(
                        make_first_order_density(D, D_a, Omega, Omega_a),
                        make_second_order_density(
                            D, D_b, D_c, D_bc, Omega, Omega_b, Omega_c, Omega_bc
                        )
                    ),
                    SymEngine::mul(
                        make_first_order_density(D, D_b, Omega, Omega_b),
                        make_second_order_density(
                            D, D_a, D_c, D_ac, Omega, Omega_a, Omega_c, Omega_ac
                        )
                    ),
                    SymEngine::mul(
                        make_first_order_density(D, D_c, Omega, Omega_c),
                        make_second_order_density(
                            D, D_a, D_b, D_ab, Omega, Omega_a, Omega_b, Omega_ab
                        )
                    )
                }))
            },
            {
                make_exc_density(D, Omega, 3),
                SymEngine::expand(SymEngine::mul({
                    make_first_order_density(D, D_a, Omega, Omega_a),
                    make_first_order_density(D, D_b, Omega, Omega_b),
                    make_first_order_density(D, D_c, Omega, Omega_c)
                }))
            }
        });
    };

    // We first do not consider derivatives of grid weights
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(
        std::string("Omega"),
        PertDependency({
            std::make_pair(a, 99),
            std::make_pair(b, 99),
            std::make_pair(c, 99),
            std::make_pair(d, 99)
        })
    );
    auto Exc_name = std::string("Exc");
    auto Exc = make_xc_energy(Exc_name, D, Omega, weight);

    REQUIRE(Exc->get_name() == Exc_name);
    REQUIRE(SymEngine::eq(*Exc->get_weight(), *weight));
    REQUIRE(SymEngine::eq(*Exc->get_state(), *D));
    REQUIRE(SymEngine::eq(*Exc->get_overlap_distribution(), *Omega));
    REQUIRE(SymEngine::unified_eq(Exc->get_weights(), SymEngine::vec_basic({weight})));
    REQUIRE(SymEngine::unified_eq(Exc->get_states(), SymEngine::vec_basic({D})));
    REQUIRE(SymEngine::unified_eq(
        Exc->get_overlap_distributions(), SymEngine::vec_basic({Omega})
    ));
    REQUIRE(Exc->get_exc_orders() == std::set<unsigned int>({0}));
    REQUIRE(eq_energy_map(
        Exc->get_energy_map(),
        ExcContractionMap({
            {weight, make_unperturbed_contraction(D, Omega)}
        })
    ));

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
    // (1) The first order XC energy density derivative
    auto Exc_a = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc->diff(a));
    REQUIRE(SymEngine::unified_eq(Exc_a->get_weights(), SymEngine::vec_basic({weight})));
    REQUIRE(SymEngine::unified_eq(Exc_a->get_states(), SymEngine::vec_basic({D, D_a})));
    REQUIRE(SymEngine::unified_eq(
        Exc_a->get_overlap_distributions(), SymEngine::vec_basic({Omega, Omega_a})
    ));
    REQUIRE(Exc_a->get_exc_orders() == std::set<unsigned int>({1}));
    REQUIRE(eq_energy_map(
        Exc_a->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                make_first_order_contraction(D, D_a, Omega, Omega_a)
            }
        })
    ));
    // (2) The second order XC energy density derivative
    auto Exc_ab = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_a->diff(b));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_weights(), SymEngine::vec_basic({weight})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_states(), SymEngine::vec_basic({D, D_b, D_a, D_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_overlap_distributions(),
        SymEngine::vec_basic({Omega, Omega_a, Omega_b, Omega_ab})
    ));
    REQUIRE(Exc_ab->get_exc_orders() == std::set<unsigned int>({1, 2}));
    REQUIRE(eq_energy_map(
        Exc_ab->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                make_second_order_contraction(
                    D, D_a, D_b, D_ab, Omega, Omega_a, Omega_b, Omega_ab
                )
            }
        })
    ));
    // (3) The third order XC energy density derivative
    auto Exc_abc = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_ab->diff(c));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_weights(), SymEngine::vec_basic({weight})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_states(),
        SymEngine::vec_basic({D, D_b, D_a, D_c, D_bc, D_ab, D_ac, D_abc})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_overlap_distributions(),
        SymEngine::vec_basic({
          Omega, Omega_c, Omega_a, Omega_b, Omega_bc, Omega_ac, Omega_ab, Omega_abc
        })
    ));
    REQUIRE(Exc_abc->get_exc_orders() == std::set<unsigned int>({1, 2, 3}));
    REQUIRE(eq_energy_map(
        Exc_abc->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                make_third_order_contraction(
                    D, D_a, D_b, D_c, D_ab, D_ac, D_bc, D_abc,
                    Omega,
                    Omega_a, Omega_b, Omega_c,
                    Omega_ab, Omega_ac, Omega_bc,
                    Omega_abc
                )
            }
        })
    ));
    // (4) The fourth order XC energy density derivative
    auto Exc_abcd = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_abc->diff(d));
    REQUIRE(SymEngine::unified_eq(
        Exc_abcd->get_weights(), SymEngine::vec_basic({weight})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abcd->get_states(),
        SymEngine::vec_basic({
          D,
          D_b, D_a, D_c, D_d,
          D_bc, D_bd, D_ab, D_ac, D_ad, D_cd,
          D_acd, D_abc, D_abd, D_bcd,
          D_abcd
        })
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abcd->get_overlap_distributions(),
        SymEngine::vec_basic({
          Omega,
          Omega_d, Omega_c, Omega_a, Omega_b,
          Omega_bd, Omega_bc, Omega_ac, Omega_ad, Omega_ab, Omega_cd,
          Omega_abd, Omega_abc, Omega_acd, Omega_bcd,
          Omega_abcd
        })
    ));
    REQUIRE(Exc_abcd->get_exc_orders() == std::set<unsigned int>({1, 2, 3, 4}));
    REQUIRE(eq_energy_map(
        Exc_abcd->get_energy_map(),
        // Equation (46)-(50), J. Chem. Phys. 140, 034103 (2014)
        ExcContractionMap({
            {
                weight,
                ExcDensityContractionMap({
                    {
                        make_exc_density(D, Omega, 1),
                        SymEngine::add({
                            make_density_vector(D, Omega_abcd),
                            make_density_vector(D_a, Omega_bcd),
                            make_density_vector(D_b, Omega_acd),
                            make_density_vector(D_c, Omega_abd),
                            make_density_vector(D_d, Omega_abc),
                            make_density_vector(D_ab, Omega_cd),
                            make_density_vector(D_ac, Omega_bd),
                            make_density_vector(D_ad, Omega_bc),
                            make_density_vector(D_bc, Omega_ad),
                            make_density_vector(D_bd, Omega_ac),
                            make_density_vector(D_cd, Omega_ab),
                            make_density_vector(D_abc, Omega_d),
                            make_density_vector(D_abd, Omega_c),
                            make_density_vector(D_acd, Omega_b),
                            make_density_vector(D_bcd, Omega_a),
                            make_density_vector(D_abcd, Omega)
                        })
                    },
                    {
                        make_exc_density(D, Omega, 2),
                        SymEngine::expand(SymEngine::add({
                            SymEngine::mul(
                                make_first_order_density(D, D_a, Omega, Omega_a),
                                make_third_order_density(
                                    D, D_b, D_c, D_d, D_bc, D_bd, D_cd, D_bcd,
                                    Omega,
                                    Omega_b, Omega_c, Omega_d,
                                    Omega_bc, Omega_bd, Omega_cd,
                                    Omega_bcd
                                )
                            ),
                            SymEngine::mul(
                                make_first_order_density(D, D_b, Omega, Omega_b),
                                make_third_order_density(
                                    D, D_a, D_c, D_d, D_ac, D_ad, D_cd, D_acd,
                                    Omega,
                                    Omega_a, Omega_c, Omega_d,
                                    Omega_ac, Omega_ad, Omega_cd,
                                    Omega_acd
                                )
                            ),
                            SymEngine::mul(
                                make_first_order_density(D, D_c, Omega, Omega_c),
                                make_third_order_density(
                                    D, D_a, D_b, D_d, D_ab, D_ad, D_bd, D_abd,
                                    Omega,
                                    Omega_a, Omega_b, Omega_d,
                                    Omega_ab, Omega_ad, Omega_bd,
                                    Omega_abd
                                )
                            ),
                            SymEngine::mul(
                                make_first_order_density(D, D_d, Omega, Omega_d),
                                make_third_order_density(
                                    D, D_a, D_b, D_c, D_ab, D_ac, D_bc, D_abc,
                                    Omega,
                                    Omega_a, Omega_b, Omega_c,
                                    Omega_ab, Omega_ac, Omega_bc,
                                    Omega_abc
                                )
                            ),
                            // The following terms are missing in equation (46),
                            // J. Chem. Phys. 140, 034103 (2014)
                            SymEngine::mul(
                                make_second_order_density(
                                    D, D_a, D_d, D_ad, Omega, Omega_a, Omega_d, Omega_ad
                                ),
                                make_second_order_density(
                                    D, D_b, D_c, D_bc, Omega, Omega_b, Omega_c, Omega_bc
                                )
                            ),
                            SymEngine::mul(
                                make_second_order_density(
                                    D, D_b, D_d, D_bd, Omega, Omega_b, Omega_d, Omega_bd
                                ),
                                make_second_order_density(
                                    D, D_a, D_c, D_ac, Omega, Omega_a, Omega_c, Omega_ac
                                )
                            ),
                            SymEngine::mul(
                                make_second_order_density(
                                    D, D_c, D_d, D_cd, Omega, Omega_c, Omega_d, Omega_cd
                                ),
                                make_second_order_density(
                                    D, D_a, D_b, D_ab, Omega, Omega_a, Omega_b, Omega_ab
                                )
                            )
                        }))
                    },
                    {
                        make_exc_density(D, Omega, 3),
                        SymEngine::expand(SymEngine::add({
                            SymEngine::mul({
                                make_first_order_density(D, D_a, Omega, Omega_a),
                                make_first_order_density(D, D_b, Omega, Omega_b),
                                make_second_order_density(
                                    D, D_c, D_d, D_cd, Omega, Omega_c, Omega_d, Omega_cd
                                )
                            }),
                            SymEngine::mul({
                                make_first_order_density(D, D_a, Omega, Omega_a),
                                make_first_order_density(D, D_c, Omega, Omega_c),
                                make_second_order_density(
                                    D, D_b, D_d, D_bd, Omega, Omega_b, Omega_d, Omega_bd
                                )
                            }),
                            SymEngine::mul({
                                make_first_order_density(D, D_a, Omega, Omega_a),
                                make_first_order_density(D, D_d, Omega, Omega_d),
                                make_second_order_density(
                                    D, D_b, D_c, D_bc, Omega, Omega_b, Omega_c, Omega_bc
                                )
                            }),
                            SymEngine::mul({
                                make_first_order_density(D, D_b, Omega, Omega_b),
                                make_first_order_density(D, D_c, Omega, Omega_c),
                                make_second_order_density(
                                    D, D_a, D_d, D_ad, Omega, Omega_a, Omega_d, Omega_ad
                                )
                            }),
                            SymEngine::mul({
                                make_first_order_density(D, D_b, Omega, Omega_b),
                                make_first_order_density(D, D_d, Omega, Omega_d),
                                make_second_order_density(
                                    D, D_a, D_c, D_ac, Omega, Omega_a, Omega_c, Omega_ac
                                )
                            }),
                            SymEngine::mul({
                                make_first_order_density(D, D_c, Omega, Omega_c),
                                make_first_order_density(D, D_d, Omega, Omega_d),
                                make_second_order_density(
                                    D, D_a, D_b, D_ab, Omega, Omega_a, Omega_b, Omega_ab
                                )
                            })
                        }))
                    },
                    {
                        make_exc_density(D, Omega, 4),
                        SymEngine::expand(SymEngine::mul({
                            make_first_order_density(D, D_a, Omega, Omega_a),
                            make_first_order_density(D, D_b, Omega, Omega_b),
                            make_first_order_density(D, D_c, Omega, Omega_c),
                            make_first_order_density(D, D_d, Omega, Omega_d)
                        }))
                    }
                })
            }
        })
    ));

    // A more complicated case is to consider the derivatives of grid weights
    weight = make_nonel_function(
        std::string("weight"),
        PertDependency({
            std::make_pair(a, 99),
            std::make_pair(b, 99)
        })
    );
    Exc = make_xc_energy(Exc_name, D, Omega, weight);
    // (1) The first order XC energy density derivative
    Exc_a = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc->diff(a));
    auto weight_a = SymEngine::rcp_dynamic_cast<const NonElecFunction>(weight->diff(a));
    REQUIRE(SymEngine::unified_eq(
        Exc_a->get_weights(), SymEngine::vec_basic({weight, weight_a})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_a->get_states(), SymEngine::vec_basic({D, D_a})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_a->get_overlap_distributions(), SymEngine::vec_basic({Omega, Omega_a})
    ));
    REQUIRE(Exc_a->get_exc_orders() == std::set<unsigned int>({0, 1}));
    REQUIRE(eq_energy_map(
        Exc_a->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                make_first_order_contraction(D, D_a, Omega, Omega_a)
            },
            {
                weight_a,
                make_unperturbed_contraction(D, Omega)
            }
        })
    ));
    // (2) The second order XC energy density derivative
    Exc_ab = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_a->diff(b));
    auto weight_b = SymEngine::rcp_dynamic_cast<const NonElecFunction>(weight->diff(b));
    auto weight_ab = SymEngine::rcp_dynamic_cast<const NonElecFunction>(weight_a->diff(b));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_weights(),
        SymEngine::vec_basic({weight, weight_b, weight_a, weight_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_states(), SymEngine::vec_basic({D, D_b, D_a, D_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_overlap_distributions(),
        SymEngine::vec_basic({Omega, Omega_a, Omega_b, Omega_ab})
    ));
    REQUIRE(Exc_ab->get_exc_orders() == std::set<unsigned int>({0, 1, 2}));
    REQUIRE(eq_energy_map(
        Exc_ab->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                make_second_order_contraction(
                    D, D_a, D_b, D_ab, Omega, Omega_a, Omega_b, Omega_ab
                )
            },
            {
                weight_a,
                make_first_order_contraction(D, D_b, Omega, Omega_b)
            },
            {
                weight_b,
                make_first_order_contraction(D, D_a, Omega, Omega_a)
            },
            {
                weight_ab,
                make_unperturbed_contraction(D, Omega)
            }
        })
    ));
    // (3) The third order XC energy density derivative, note that the grid
    // weights do not depend on the perturbation `c`
    Exc_abc = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_ab->diff(c));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_weights(),
        SymEngine::vec_basic({weight, weight_b, weight_a, weight_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_states(),
        SymEngine::vec_basic({D, D_b, D_a, D_c, D_bc, D_ab, D_ac, D_abc})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_overlap_distributions(),
        SymEngine::vec_basic({
          Omega, Omega_c, Omega_a, Omega_b, Omega_bc, Omega_ac, Omega_ab, Omega_abc
        })
    ));
    REQUIRE(Exc_abc->get_exc_orders() == std::set<unsigned int>({1, 2, 3}));
    REQUIRE(eq_energy_map(
        Exc_abc->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                // Exactly equation (45), J. Chem. Phys. 140, 034103 (2014)
                make_third_order_contraction(
                    D,
                    D_a, D_b, D_c,
                    D_ab, D_ac, D_bc,
                    D_abc,
                    Omega,
                    Omega_a, Omega_b, Omega_c,
                    Omega_ab, Omega_ac, Omega_bc,
                    Omega_abc
                )
            },
            {
                weight_a,
                make_second_order_contraction(
                    D, D_b, D_c, D_bc, Omega, Omega_b, Omega_c, Omega_bc
                )
            },
            {
                weight_b,
                make_second_order_contraction(
                    D, D_a, D_c, D_ac, Omega, Omega_a, Omega_c, Omega_ac
                )
            },
            {
                weight_ab,
                make_first_order_contraction(D, D_c, Omega, Omega_c)
            }
        })
    ));

    // Last example is when the generalized overlap distributions do not depend
    // on all perturbations
    Omega = make_1el_operator(
        std::string("Omega"),
        PertDependency({
            std::make_pair(c, 99),
            std::make_pair(d, 99)
        })
    );
    Omega_c = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(c));
    Omega_d = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(d));
    Omega_cd = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_c->diff(d));
    Exc = make_xc_energy(Exc_name, D, Omega, weight);
    // (1) The first order XC energy density derivative
    Exc_a = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc->diff(a));
    REQUIRE(SymEngine::unified_eq(
        Exc_a->get_weights(), SymEngine::vec_basic({weight, weight_a})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_a->get_states(), SymEngine::vec_basic({D, D_a})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_a->get_overlap_distributions(), SymEngine::vec_basic({Omega})
    ));
    REQUIRE(Exc_a->get_exc_orders() == std::set<unsigned int>({0, 1}));
    REQUIRE(eq_energy_map(
        Exc_a->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                ExcDensityContractionMap({
                    {make_exc_density(D, Omega, 1), make_density_vector(D_a, Omega)}
                })
            },
            {
                weight_a,
                make_unperturbed_contraction(D, Omega)
            }
        })
    ));
    // (2) The second order XC energy density derivative
    Exc_ab = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_a->diff(b));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_weights(),
        SymEngine::vec_basic({weight, weight_b, weight_a, weight_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_states(), SymEngine::vec_basic({D, D_b, D_a, D_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_ab->get_overlap_distributions(), SymEngine::vec_basic({Omega})
    ));
    REQUIRE(Exc_ab->get_exc_orders() == std::set<unsigned int>({0, 1, 2}));
    REQUIRE(eq_energy_map(
        Exc_ab->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                ExcDensityContractionMap({
                    {
                        make_exc_density(D, Omega, 1),
                        make_density_vector(D_ab, Omega)
                    },
                    {
                        make_exc_density(D, Omega, 2),
                        SymEngine::expand(SymEngine::mul(
                            make_density_vector(D_a, Omega),
                            make_density_vector(D_b, Omega)
                        ))
                    }
                })
            },
            {
                weight_a,
                ExcDensityContractionMap({
                    {make_exc_density(D, Omega, 1), make_density_vector(D_b, Omega)}
                })
            },
            {
                weight_b,
                ExcDensityContractionMap({
                    {make_exc_density(D, Omega, 1), make_density_vector(D_a, Omega)}
                })
            },
            {
                weight_ab,
                make_unperturbed_contraction(D, Omega)
            }
        })
    ));
    // (3) The third order XC energy density derivative, note that the grid
    // weights do not depend on the perturbation `c`
    Exc_abc = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_ab->diff(c));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_weights(),
        SymEngine::vec_basic({weight, weight_b, weight_a, weight_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_states(),
        SymEngine::vec_basic({D, D_b, D_a, D_c, D_bc, D_ab, D_ac, D_abc})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abc->get_overlap_distributions(),
        SymEngine::vec_basic({Omega, Omega_c})
    ));
    REQUIRE(Exc_abc->get_exc_orders() == std::set<unsigned int>({1, 2, 3}));
    REQUIRE(eq_energy_map(
        Exc_abc->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                // By removing derivatives of `Omega` with respect to `a` and
                // `b` from `make_third_order_contraction()`
                ExcDensityContractionMap({
                    {
                        make_exc_density(D, Omega, 1),
                        make_first_order_density(D_ab, D_abc, Omega, Omega_c)
                    },
                    {
                        make_exc_density(D, Omega, 2),
                        SymEngine::expand(SymEngine::add({
                            SymEngine::mul(
                                make_density_vector(D_a, Omega),
                                make_first_order_density(D_b, D_bc, Omega, Omega_c)
                            ),
                            SymEngine::mul(
                                make_density_vector(D_b, Omega),
                                make_first_order_density(D_a, D_ac, Omega, Omega_c)
                            ),
                            SymEngine::mul(
                                make_first_order_density(D, D_c, Omega, Omega_c),
                                make_density_vector(D_ab, Omega)
                            )
                        }))
                    },
                    {
                        make_exc_density(D, Omega, 3),
                        SymEngine::expand(SymEngine::mul({
                            make_density_vector(D_a, Omega),
                            make_density_vector(D_b, Omega),
                            make_first_order_density(D, D_c, Omega, Omega_c)
                        }))
                    }
                })
            },
            {
                weight_a,
                ExcDensityContractionMap({
                    {
                        make_exc_density(D, Omega, 1),
                        make_first_order_density(D_b, D_bc, Omega, Omega_c)
                    },
                    {
                        make_exc_density(D, Omega, 2),
                        SymEngine::expand(SymEngine::mul(
                            make_density_vector(D_b, Omega),
                            make_first_order_density(D, D_c, Omega, Omega_c)
                        ))
                    }
                })
            },
            {
                weight_b,
                ExcDensityContractionMap({
                    {
                        make_exc_density(D, Omega, 1),
                        make_first_order_density(D_a, D_ac, Omega, Omega_c)
                    },
                    {
                        make_exc_density(D, Omega, 2),
                        SymEngine::expand(SymEngine::mul(
                            make_density_vector(D_a, Omega),
                            make_first_order_density(D, D_c, Omega, Omega_c)
                        ))
                    }
                })
            },
            {
                weight_ab,
                make_first_order_contraction(D, D_c, Omega, Omega_c)
            }
        })
    ));
    // (4) The fourth order XC energy density derivative
    Exc_abcd = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_abc->diff(d));
    REQUIRE(SymEngine::unified_eq(
        Exc_abcd->get_weights(),
        SymEngine::vec_basic({weight, weight_b, weight_a, weight_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abcd->get_states(),
        SymEngine::vec_basic({
          D,
          D_b, D_a, D_c, D_d,
          D_bc, D_bd, D_ab, D_ac, D_ad, D_cd,
          D_acd, D_abc, D_abd, D_bcd,
          D_abcd
        })
    ));
    REQUIRE(SymEngine::unified_eq(
        Exc_abcd->get_overlap_distributions(),
        SymEngine::vec_basic({Omega, Omega_d, Omega_c, Omega_cd})
    ));
    REQUIRE(Exc_abcd->get_exc_orders() == std::set<unsigned int>({1, 2, 3, 4}));
    REQUIRE(eq_energy_map(
        Exc_abcd->get_energy_map(),
        ExcContractionMap({
            {
                weight,
                ExcDensityContractionMap({
                    {
                        make_exc_density(D, Omega, 1),
                        make_second_order_density(
                            D_ab, D_abc, D_abd, D_abcd, Omega, Omega_c, Omega_d, Omega_cd
                        )
                    },
                    {
                        make_exc_density(D, Omega, 2),
                        SymEngine::expand(SymEngine::add({
                            SymEngine::mul(
                                make_first_order_density(D_ab, D_abc, Omega, Omega_c),
                                make_first_order_density(D, D_d, Omega, Omega_d)
                            ),
                            SymEngine::mul(
                                make_first_order_density(D_a, D_ad, Omega, Omega_d),
                                make_first_order_density(D_b, D_bc, Omega, Omega_c)
                            ),
                            SymEngine::mul(
                                make_density_vector(D_a, Omega),
                                make_second_order_density(
                                    D_b, D_bc, D_bd, D_bcd,
                                    Omega, Omega_c, Omega_d, Omega_cd
                                )
                            ),
                            SymEngine::mul(
                                make_first_order_density(D_b, D_bd, Omega, Omega_d),
                                make_first_order_density(D_a, D_ac, Omega, Omega_c)
                            ),
                            SymEngine::mul(
                                make_density_vector(D_b, Omega),
                                make_second_order_density(
                                    D_a, D_ac, D_ad, D_acd,
                                    Omega, Omega_c, Omega_d, Omega_cd
                                )
                            ),
                            SymEngine::mul(
                                make_density_vector(D_ab, Omega),
                                make_second_order_density(
                                    D, D_c, D_d, D_cd, Omega, Omega_c, Omega_d, Omega_cd
                                )
                            ),
                            SymEngine::mul(
                                make_first_order_density(D, D_c, Omega, Omega_c),
                                make_first_order_density(D_ab, D_abd, Omega, Omega_d)
                            )
                        }))
                    },
                    {
                        make_exc_density(D, Omega, 3),
                        SymEngine::expand(SymEngine::add({
                            SymEngine::mul({
                                make_density_vector(D_a, Omega),
                                make_first_order_density(D_b, D_bc, Omega, Omega_c),
                                make_first_order_density(D, D_d, Omega, Omega_d)
                            }),
                            SymEngine::mul({
                                make_density_vector(D_b, Omega),
                                make_first_order_density(D_a, D_ac, Omega, Omega_c),
                                make_first_order_density(D, D_d, Omega, Omega_d)
                            }),
                            SymEngine::mul({
                                make_density_vector(D_ab, Omega),
                                make_first_order_density(D, D_c, Omega, Omega_c),
                                make_first_order_density(D, D_d, Omega, Omega_d)
                            }),
                            SymEngine::mul({
                                make_density_vector(D_b, Omega),
                                make_first_order_density(D_a, D_ad, Omega, Omega_d),
                                make_first_order_density(D, D_c, Omega, Omega_c)
                            }),
                            SymEngine::mul({
                                make_density_vector(D_a, Omega),
                                make_first_order_density(D_b, D_bd, Omega, Omega_d),
                                make_first_order_density(D, D_c, Omega, Omega_c)
                            }),
                            SymEngine::mul({
                                make_density_vector(D_a, Omega),
                                make_density_vector(D_b, Omega),
                                make_second_order_density(
                                    D, D_c, D_d, D_cd, Omega, Omega_c, Omega_d, Omega_cd
                                )
                            })
                        }))
                    },
                    {
                        make_exc_density(D, Omega, 4),
                        SymEngine::expand(SymEngine::mul({
                            make_density_vector(D_a, Omega),
                            make_density_vector(D_b, Omega),
                            make_first_order_density(D, D_c, Omega, Omega_c),
                            make_first_order_density(D, D_d, Omega, Omega_d),
                        }))
                    }
                })
            },
            {
                weight_a,
                ExcDensityContractionMap({
                    {
                        make_exc_density(D, Omega, 1),
                        make_second_order_density(
                            D_b, D_bc, D_bd, D_bcd, Omega, Omega_c, Omega_d, Omega_cd
                        )
                    },
                    {
                        make_exc_density(D, Omega, 2),
                        SymEngine::expand(SymEngine::add({
                            SymEngine::mul(
                                make_first_order_density(D, D_d, Omega, Omega_d),
                                make_first_order_density(D_b, D_bc, Omega, Omega_c)
                            ),
                            SymEngine::mul(
                                make_first_order_density(D_b, D_bd, Omega, Omega_d),
                                make_first_order_density(D, D_c, Omega, Omega_c)
                            ),
                            SymEngine::mul(
                                make_density_vector(D_b, Omega),
                                make_second_order_density(
                                    D, D_c, D_d, D_cd, Omega, Omega_c, Omega_d, Omega_cd
                                )
                            )
                        }))
                    },
                    {
                        make_exc_density(D, Omega, 3),
                        SymEngine::expand(SymEngine::mul({
                            make_density_vector(D_b, Omega),
                            make_first_order_density(D, D_c, Omega, Omega_c),
                            make_first_order_density(D, D_d, Omega, Omega_d)
                        }))
                    }
                })
            },
            {
                weight_b,
                ExcDensityContractionMap({
                    {
                        make_exc_density(D, Omega, 1),
                        make_second_order_density(
                            D_a, D_ac, D_ad, D_acd, Omega, Omega_c, Omega_d, Omega_cd
                        )
                    },
                    {
                        make_exc_density(D, Omega, 2),
                        SymEngine::expand(SymEngine::add({
                            SymEngine::mul(
                                make_first_order_density(D, D_d, Omega, Omega_d),
                                make_first_order_density(D_a, D_ac, Omega, Omega_c)
                            ),
                            SymEngine::mul(
                                make_first_order_density(D_a, D_ad, Omega, Omega_d),
                                make_first_order_density(D, D_c, Omega, Omega_c)
                            ),
                            SymEngine::mul(
                                make_density_vector(D_a, Omega),
                                make_second_order_density(
                                    D, D_c, D_d, D_cd, Omega, Omega_c, Omega_d, Omega_cd
                                )
                            )
                        }))
                    },
                    {
                        make_exc_density(D, Omega, 3),
                        SymEngine::expand(SymEngine::mul({
                            make_density_vector(D_a, Omega),
                            make_first_order_density(D, D_c, Omega, Omega_c),
                            make_first_order_density(D, D_d, Omega, Omega_d)
                        }))
                    }
                })
            },
            {
                weight_ab,
                make_second_order_contraction(
                    D, D_c, D_d, D_cd, Omega, Omega_c, Omega_d, Omega_cd
                )
            }
        })
    ));
}

TEST_CASE("Test ExchCorrPotential and make_xc_potential()", "[ExchCorrPotential]")
{
    auto D = make_1el_density(std::string("D"));
    auto a = make_perturbation(std::string("a"));
    auto b = make_perturbation(std::string("b"));

    // Lambda expressions to make contractions of (un)perturbed generalized
    // density vectors and the XC functional derivative vectors
    //
    // For unperturbed XC potential, there is only one unperturbed grid
    // weight, and order of the XC energy functional derivative is 1, and no
    // generalized density vectors to be contracted with
    auto make_unperturbed_contraction = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const OneElecOperator>& Omega
    ) {
        return ExcDensityContractionMap({
            {make_exc_density(D, Omega, 1), SymEngine::RCP<const SymEngine::Basic>()}
        });
    };
    // Equation (47), J. Chem. Phys. 140, 034103 (2014)
    auto make_first_order_density = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a
    ) {
        return SymEngine::add(
            make_density_vector(D, Omega_a),
            make_density_vector(D_a, Omega)
        );
    };
    // Equation (59), J. Chem. Phys. 140, 034103 (2014)
    auto make_first_order_contraction = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a
    ) {
        return ExcDensityContractionMap({
            {
                make_exc_density(D, Omega, 2),
                make_first_order_density(D, D_a, Omega, Omega_a)
            }
        });
    };
    // Equation (48), J. Chem. Phys. 140, 034103 (2014)
    auto make_second_order_density = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const ElectronicState>& D_b,
        const SymEngine::RCP<const ElectronicState>& D_ab,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a,
        const SymEngine::RCP<const OneElecOperator>& Omega_b,
        const SymEngine::RCP<const OneElecOperator>& Omega_ab
    ) {
        return SymEngine::add({
            make_density_vector(D, Omega_ab),
            make_density_vector(D_a, Omega_b),
            make_density_vector(D_b, Omega_a),
            make_density_vector(D_ab, Omega)
        });
    };
    // Equation (60), J. Chem. Phys. 140, 034103 (2014)
    auto make_second_order_contraction = [&](
        const SymEngine::RCP<const ElectronicState>& D,
        const SymEngine::RCP<const ElectronicState>& D_a,
        const SymEngine::RCP<const ElectronicState>& D_b,
        const SymEngine::RCP<const ElectronicState>& D_ab,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const OneElecOperator>& Omega_a,
        const SymEngine::RCP<const OneElecOperator>& Omega_b,
        const SymEngine::RCP<const OneElecOperator>& Omega_ab
    ) {
        return ExcDensityContractionMap({
            {
                make_exc_density(D, Omega, 2),
                make_second_order_density(
                    D, D_a, D_b, D_ab, Omega, Omega_a, Omega_b, Omega_ab
                )
            },
            {
                make_exc_density(D, Omega, 3),
                SymEngine::expand(SymEngine::mul(
                    make_first_order_density(D, D_a, Omega, Omega_a),
                    make_first_order_density(D, D_b, Omega, Omega_b)
                ))
            }
        });
    };

    // We first do not consider derivatives of grid weights
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(
        std::string("Omega"),
        PertDependency({std::make_pair(a, 99), std::make_pair(b, 99)})
    );
    auto Vxc_name = std::string("Vxc");
    auto Vxc = make_xc_potential(Vxc_name, D, Omega, weight);

    REQUIRE(Vxc->get_name() == Vxc_name);
    REQUIRE(SymEngine::eq(*Vxc->get_weight(), *weight));
    REQUIRE(SymEngine::eq(*Vxc->get_state(), *D));
    REQUIRE(SymEngine::eq(*Vxc->get_overlap_distribution(), *Omega));
    REQUIRE(SymEngine::unified_eq(
        Vxc->get_weights(), SymEngine::vec_basic({weight})
    ));
    REQUIRE(SymEngine::unified_eq(Vxc->get_states(), SymEngine::vec_basic({D})));
    REQUIRE(SymEngine::unified_eq(
        Vxc->get_overlap_distributions(), SymEngine::vec_basic({Omega})
    ));
    REQUIRE(Vxc->get_exc_orders() == std::set<unsigned int>({1}));
    REQUIRE(eq_potential_map(
        Vxc->get_potential_map(),
        VxcContractionMap({
            {
                Omega,
                ExcContractionMap({{weight, make_unperturbed_contraction(D, Omega)}})
            }
        })
    ));

    // Tests from J. Chem. Phys. 140, 034103 (2014), equations (57)-(60).
    auto D_a = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(a));
    auto D_b = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(b));
    auto D_ab = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_a->diff(b));
    auto Omega_a = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(a));
    auto Omega_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(b));
    auto Omega_ab = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_a->diff(b));
    // (1) The first order XC potential matrix
    auto Vxc_a = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Vxc->diff(a));
    REQUIRE(SymEngine::unified_eq(
        Vxc_a->get_weights(), SymEngine::vec_basic({weight})
    ));
    REQUIRE(SymEngine::unified_eq(
        Vxc_a->get_states(), SymEngine::vec_basic({D, D_a})
    ));
    REQUIRE(SymEngine::unified_eq(
        Vxc_a->get_overlap_distributions(), SymEngine::vec_basic({Omega, Omega_a})
    ));
    REQUIRE(Vxc_a->get_exc_orders() == std::set<unsigned int>({1, 2}));
    REQUIRE(eq_potential_map(
        Vxc_a->get_potential_map(),
        VxcContractionMap({
            {
                Omega,
                ExcContractionMap({
                    {weight, make_first_order_contraction(D, D_a, Omega, Omega_a)}
                })
            },
            {
                Omega_a,
                ExcContractionMap({{weight, make_unperturbed_contraction(D, Omega)}})
            }
        })
    ));
    // (2) The first order XC potential matrix
    auto Vxc_ab = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Vxc_a->diff(b));
    REQUIRE(SymEngine::unified_eq(
        Vxc_ab->get_weights(), SymEngine::vec_basic({weight})
    ));
    REQUIRE(SymEngine::unified_eq(
        Vxc_ab->get_states(), SymEngine::vec_basic({D, D_b, D_a, D_ab})
    ));
    REQUIRE(SymEngine::unified_eq(
        Vxc_ab->get_overlap_distributions(),
        SymEngine::vec_basic({Omega, Omega_a, Omega_b, Omega_ab})
    ));
    REQUIRE(Vxc_ab->get_exc_orders() == std::set<unsigned int>({1, 2, 3}));
    REQUIRE(eq_potential_map(
        Vxc_ab->get_potential_map(),
        VxcContractionMap({
            {
                Omega,
                ExcContractionMap({
                    {
                        weight,
                        make_second_order_contraction(
                            D, D_a, D_b, D_ab, Omega, Omega_a, Omega_b, Omega_ab
                        )
                    }
                })
            },
            {
                Omega_a,
                ExcContractionMap({
                    {weight, make_first_order_contraction(D, D_b, Omega, Omega_b)}
                })
            },
            {
                Omega_b,
                ExcContractionMap({
                    {weight, make_first_order_contraction(D, D_a, Omega, Omega_a)}
                })
            },
            {
                Omega_ab,
                ExcContractionMap({{weight, make_unperturbed_contraction(D, Omega)}})
            }
        })
    ));

    // A more complicated case is to consider the derivatives of grid weights

    // Last example is when the generalized overlap distributions do not depend
    // on all perturbations
}
