#define CATCH_CONFIG_MAIN

#include <string>

#include <catch2/catch.hpp>

#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"

using namespace Tinned;

// Equation (80), J. Chem. Phys. 129, 214108 (2008)
inline SymEngine::RCP<const SymEngine::Basic> make_ks_energy(
    const SymEngine::RCP<const OneElecOperator>& h,
    const SymEngine::RCP<const OneElecOperator>& V,
    const SymEngine::RCP<const TwoElecOperator>& G,
    const SymEngine::RCP<const OneElecDensity>& D,
    const SymEngine::RCP<const ExchCorrEnergy>& Exc,
    const SymEngine::RCP<const NonElecFunction>& hnuc
)
{
    return SymEngine::add({
        SymEngine::trace(SymEngine::matrix_mul({h, D})),
        SymEngine::trace(SymEngine::matrix_mul({V, D})),
        make_2el_energy(G),
        Exc,
        hnuc
    });
}

// Equation (229), J. Chem. Phys. 129, 214108 (2008)
inline SymEngine::RCP<const SymEngine::Basic> make_tdscf_equation(
    const SymEngine::RCP<const SymEngine::Basic>& F,
    const SymEngine::RCP<const OneElecDensity>& D,
    const SymEngine::RCP<const OneElecOperator>& S
)
{
    auto Dt = make_dt_operator(D);
    auto St = make_dt_operator(S);
    auto minus_one_half = SymEngine::div(SymEngine::minus_one, SymEngine::two);
    return SymEngine::matrix_add({
        SymEngine::matrix_mul({F, D, S}),
        SymEngine::matrix_mul({SymEngine::minus_one, S, D, F}),
        SymEngine::matrix_mul({SymEngine::minus_one, S, Dt, S}),
        SymEngine::matrix_mul({minus_one_half, St, D, S}),
        SymEngine::matrix_mul({minus_one_half, S, D, St})
    });
}

TEST_CASE("Test KeepVisitor and keep_if()", "[KeepVisitor]")
{
    auto a = make_perturbation(std::string("a"));
    auto b = make_perturbation(std::string("b"));
    auto dependencies = PertDependency({std::make_pair(a, 99), std::make_pair(b, 99)});
    auto D = make_1el_density(std::string("D"));
    auto h = make_1el_operator(std::string("h"), dependencies);
    auto V = make_1el_operator(std::string("V"), dependencies);
    auto G = make_2el_operator(std::string("G"), D, dependencies);
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(std::string("Omega"), dependencies);
    auto Exc = make_xc_energy(std::string("Exc"), D, Omega, weight);
    auto Fxc = make_xc_potential(std::string("Fxc"), D, Omega, weight);
    auto hnuc = make_nonel_function(std::string("hnuc"), dependencies);
    // Equation (80), J. Chem. Phys. 129, 214108 (2008)
    auto E = make_ks_energy(h, V, G, D, Exc, hnuc);
    auto F = SymEngine::matrix_add({h, G, V, Fxc});

    REQUIRE(SymEngine::eq(
        *keep_if(E, SymEngine::set_basic({h})),
        *SymEngine::trace(SymEngine::matrix_mul({h, D}))
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(E, SymEngine::set_basic({h, V})),
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({h, D})),
            SymEngine::trace(SymEngine::matrix_mul({V, D}))
        })
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(E, SymEngine::set_basic({D})),
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({h, D})),
            SymEngine::trace(SymEngine::matrix_mul({V, D})),
            make_2el_energy(G),
            Exc
        })
    ));

    // The first order
    auto E_a = E->diff(a);
    auto F_a = F->diff(a);
    auto D_a = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(a));
    auto h_a = SymEngine::rcp_dynamic_cast<const OneElecOperator>(h->diff(a));
    auto V_a = SymEngine::rcp_dynamic_cast<const OneElecOperator>(V->diff(a));
    //FIXME: The following construction of `Ga` is not recommended because
    // `derivatives` are not checked against `dependencies`
    auto Ga = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({a})
    );
    auto G_Da = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_a, dependencies
    );
    auto Exc_a = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc->diff(a));
    auto Exc_Da = SymEngine::make_rcp<const ExchCorrEnergy>(
        Exc->get_name(),
        Exc->get_state(),
        Exc->get_overlap_distribution(),
        Exc->get_weight(),
        SymEngine::mul({
            weight,
            make_exc_density(D, Omega, 1),
            SymEngine::trace(SymEngine::matrix_mul({Omega, D_a}))
        })
    );
    auto hnuc_a = SymEngine::rcp_dynamic_cast<const NonElecFunction>(hnuc->diff(a));
    auto Fxc_a = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Fxc->diff(a));
    REQUIRE(keep_if(E_a, SymEngine::set_basic({hnuc})).is_null());
    REQUIRE(SymEngine::eq(
        *keep_if(E_a, SymEngine::set_basic({h})),
        *SymEngine::trace(SymEngine::matrix_mul({h, D_a}))
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(E_a, SymEngine::set_basic({h_a})),
        *SymEngine::trace(SymEngine::matrix_mul({h_a, D}))
    ));
    // Unperturbed two-electron operator contracted with perturbed density matrix
    auto E2_Da = make_2el_energy(G, D_a);
    REQUIRE(SymEngine::eq(
        *keep_if(E_a, SymEngine::set_basic({G})),
        *SymEngine::mul(SymEngine::two, E2_Da)
    ));
    // Perturbed two-electron operator contracted with unperturbed density matrix
    REQUIRE(SymEngine::eq(
        *keep_if(E_a, SymEngine::set_basic({Ga})),
        *make_2el_energy(Ga, D)
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(E_a, SymEngine::set_basic({D})),
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({h_a, D})),
            SymEngine::trace(SymEngine::matrix_mul({V_a, D})),
            make_2el_energy(G)->diff(a),
            Exc_a
        })
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(E_a, SymEngine::set_basic({D_a})),
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({h, D_a})),
            SymEngine::trace(SymEngine::matrix_mul({V, D_a})),
            E2_Da,
            make_2el_energy(G_Da, D),
            Exc_Da
        })
    ));
    REQUIRE(SymEngine::eq(
        *E_a,
        *SymEngine::add(keep_if(E_a, SymEngine::set_basic({D, D_a})), hnuc_a)
    ));
    REQUIRE(keep_if(F_a, SymEngine::set_basic({h})).is_null());
    REQUIRE(SymEngine::eq(*keep_if(F_a, SymEngine::set_basic({h_a})), *h_a));
    REQUIRE(SymEngine::eq(
        *keep_if(F_a, SymEngine::set_basic({D})),
        *SymEngine::matrix_add({Ga, Fxc_a})
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(F_a, SymEngine::set_basic({D_a})),
        *SymEngine::matrix_add({
            G_Da,
            SymEngine::make_rcp<const ExchCorrPotential>(
                Fxc->get_name(),
                Fxc->get_state(),
                Fxc->get_overlap_distribution(),
                Fxc->get_weight(),
                // The first term of Equation (57), J. Chem. Phys. 140, 034103
                // (2014), but only keeping the second term of Equation (47)
                SymEngine::matrix_mul({
                    SymEngine::mul({
                        weight,
                        make_exc_density(D, Omega, 2),
                        make_density_vector(D_a, Omega)
                    }),
                    Omega
                })
            )
        })
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(F_a, SymEngine::set_basic({D, D_a})),
        *SymEngine::matrix_add({G->diff(a), Fxc_a})
    ));

    // The second order
    auto E_ab = E_a->diff(b);
    auto F_ab = F_a->diff(b);
    auto D_b = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(b));
    auto D_ab = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_a->diff(b));
    auto h_ab = SymEngine::rcp_dynamic_cast<const OneElecOperator>(h_a->diff(b));
    auto V_ab = SymEngine::rcp_dynamic_cast<const OneElecOperator>(V_a->diff(b));
    auto Gb = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({b})
    );
    auto Gb_Da = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_a, dependencies, SymEngine::multiset_basic({b})
    );
    auto G_Dab = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_ab, dependencies
    );
    auto Ga_Db = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_b, dependencies, SymEngine::multiset_basic({a})
    );
    auto Gab = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({a, b})
    );
    auto Omega_a = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(a));
    auto Omega_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(b));
    auto Exc_ab = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_a->diff(b));
    auto Fxc_ab = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Fxc_a->diff(b));
    REQUIRE(SymEngine::eq(
        *keep_if(E_ab, SymEngine::set_basic({h_a})),
        *SymEngine::trace(SymEngine::matrix_mul({h_a, D_b}))
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(E_ab, SymEngine::set_basic({D})),
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({h_ab, D})),
            SymEngine::trace(SymEngine::matrix_mul({V_ab, D})),
            make_2el_energy(Gb, D_a),
            make_2el_energy(G, D_ab),
            make_2el_energy(Gb_Da, D),
            make_2el_energy(G_Dab, D),
            make_2el_energy(Ga_Db, D),
            make_2el_energy(Ga, D_b),
            make_2el_energy(Gab, D),
            Exc_ab
        })
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(F_ab, SymEngine::set_basic({D})), *SymEngine::matrix_add({Gab, Fxc_ab})
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(F_ab, SymEngine::set_basic({D_a})),
        *SymEngine::matrix_add({
            Gb_Da,
            SymEngine::make_rcp<const ExchCorrPotential>(
                Fxc->get_name(),
                Fxc->get_state(),
                Fxc->get_overlap_distribution(),
                Fxc->get_weight(),
                SymEngine::matrix_add({
                    // The first term of Equation (58), J. Chem. Phys. 140,
                    // 034103 (2014), but only keeping the second term of
                    // Equation (47) and the third term of Equation (48)
                    SymEngine::matrix_mul({
                        SymEngine::add(
                            SymEngine::mul({
                                weight,
                                make_exc_density(D, Omega, 3),
                                make_density_vector(D_a, Omega),
                                (make_density_vector(D, Omega))->diff(b)
                            }),
                            SymEngine::mul({
                                weight,
                                make_exc_density(D, Omega, 2),
                                make_density_vector(D_a, Omega_b)
                            })
                        ),
                        Omega
                    }),
                    // The second term of Equation (58), J. Chem. Phys. 140,
                    // 034103 (2014), but only keeping the second term of
                    // Equation (47)
                    SymEngine::matrix_mul({
                        SymEngine::mul({
                            weight,
                            make_exc_density(D, Omega, 2),
                            make_density_vector(D_a, Omega)
                        }),
                        Omega_b
                    })
                })
            )
        })
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(F_ab, SymEngine::set_basic({D_ab})),
        *SymEngine::matrix_add({
            G_Dab,
            SymEngine::make_rcp<const ExchCorrPotential>(
                Fxc->get_name(),
                Fxc->get_state(),
                Fxc->get_overlap_distribution(),
                Fxc->get_weight(),
                // The first term of Equation (58), J. Chem. Phys. 140, 034103
                // (2014), but only keeping the last term of Equation (48)
                SymEngine::matrix_mul({
                    SymEngine::mul({
                        weight,
                        make_exc_density(D, Omega, 2),
                        make_density_vector(D_ab, Omega)
                    }),
                    Omega
                })
            )
        })
    ));
    REQUIRE(SymEngine::eq(
        *keep_if(F_ab, SymEngine::set_basic({D_a, D_b})),
        *SymEngine::matrix_add({
            Ga_Db,
            Gb_Da,
            SymEngine::make_rcp<const ExchCorrPotential>(
                Fxc->get_name(),
                Fxc->get_state(),
                Fxc->get_overlap_distribution(),
                Fxc->get_weight(),
                SymEngine::matrix_add({
                    // The first term of Equation (58), J. Chem. Phys. 140,
                    // 034103 (2014), but removing the product of unperturbed
                    // density matrix of Equation (47) and keeping the second
                    // and the third terms of Equation (48)
                    SymEngine::matrix_mul({
                        SymEngine::add(
                            SymEngine::mul({
                                weight,
                                make_exc_density(D, Omega, 3),
                                SymEngine::add({
                                    SymEngine::mul(
                                        make_density_vector(D, Omega_a),
                                        make_density_vector(D_b, Omega)
                                    ),
                                    SymEngine::mul(
                                        make_density_vector(D_a, Omega),
                                        make_density_vector(D, Omega_b)
                                    ),
                                    SymEngine::mul(
                                        make_density_vector(D_a, Omega),
                                        make_density_vector(D_b, Omega)
                                    )
                                })
                            }),
                            SymEngine::mul({
                                weight,
                                make_exc_density(D, Omega, 2),
                                SymEngine::add(
                                    make_density_vector(D_a, Omega_b),
                                    make_density_vector(D_b, Omega_a)
                                )
                            })
                        ),
                        Omega
                    }),
                    // The second term of Equation (58), J. Chem. Phys. 140,
                    // 034103 (2014), but only keeping the second term of
                    // Equation (47)
                    SymEngine::matrix_mul({
                        SymEngine::mul({
                            weight,
                            make_exc_density(D, Omega, 2),
                            make_density_vector(D_a, Omega)
                        }),
                        Omega_b
                    }),
                    // The third term of Equation (58), J. Chem. Phys. 140,
                    // 034103 (2014), but only keeping the second term of
                    // Equation (47)
                    SymEngine::matrix_mul({
                        SymEngine::mul({
                            weight,
                            make_exc_density(D, Omega, 2),
                            make_density_vector(D_b, Omega)
                        }),
                        Omega_a
                    })
                })
            )
        })
    ));
}

TEST_CASE("Test RemoveVisitor and remove_if()", "[RemoveVisitor]")
{
    auto a = make_perturbation(std::string("a"));
    auto b = make_perturbation(std::string("b"));
    auto dependencies = PertDependency({std::make_pair(a, 99), std::make_pair(b, 99)});
    auto D = make_1el_density(std::string("D"));
    auto h = make_1el_operator(std::string("h"), dependencies);
    auto V = make_1el_operator(std::string("V"), dependencies);
    auto G = make_2el_operator(std::string("G"), D, dependencies);
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(std::string("Omega"), dependencies);
    auto Exc = make_xc_energy(std::string("Exc"), D, Omega, weight);
    auto hnuc = make_nonel_function(std::string("hnuc"), dependencies);
    // Equation (80), J. Chem. Phys. 129, 214108 (2008)
    auto E = make_ks_energy(h, V, G, D, Exc, hnuc);
    auto E_a = E->diff(a);
    auto D_a = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(a));
    auto h_a = SymEngine::rcp_dynamic_cast<const OneElecOperator>(h->diff(a));
    auto V_a = SymEngine::rcp_dynamic_cast<const OneElecOperator>(V->diff(a));
    auto Ga = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({a})
    );
    // weight*exc^(1)(tr(Omega*D))*tr(Derivative(Omega, a)*D)
    // Equations (43) and (47) by removing `D_a`, J. Chem. Phys. 140, 034103 (2014)
    auto Omega_a = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(a));
    auto Exc_a_D = SymEngine::make_rcp<const ExchCorrEnergy>(
        Exc->get_name(),
        Exc->get_state(),
        Exc->get_overlap_distribution(),
        Exc->get_weight(),
        SymEngine::mul({
            weight,
            make_exc_density(D, Omega, 1),
            SymEngine::trace(SymEngine::matrix_mul({Omega_a, D}))
        })
    );
    auto hnuc_a = SymEngine::rcp_dynamic_cast<const NonElecFunction>(hnuc->diff(a));
    // Equation (81), J. Chem. Phys. 129, 214108 (2008)
    auto E_0a = remove_if(E_a, SymEngine::set_basic({D_a}));
    REQUIRE(SymEngine::eq(
        *E_0a, *make_ks_energy(h_a, V_a, Ga, D, Exc_a_D, hnuc_a)
    ));
    auto E_ab = E_a->diff(b);
    auto D_b = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(b));
    auto D_ab = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_a->diff(b));
    auto h_ab = SymEngine::rcp_dynamic_cast<const OneElecOperator>(h_a->diff(b));
    auto V_ab = SymEngine::rcp_dynamic_cast<const OneElecOperator>(V_a->diff(b));
    auto Ga_Db = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_b, dependencies, SymEngine::multiset_basic({a})
    );
    auto Gab = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({a, b})
    );
    // Equations (44), (47) and (48) by removing `D_a` and `D_ab`,
    // J. Chem. Phys. 140, 034103 (2014), which contains all XC contributions
    // from Equations (208) and (209), J. Chem. Phys. 129, 214108 (2008)
    auto Omega_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(b));
    auto Omega_ab = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_a->diff(b));
    auto Exc_ab_D = SymEngine::make_rcp<const ExchCorrEnergy>(
        Exc->get_name(),
        Exc->get_state(),
        Exc->get_overlap_distribution(),
        Exc->get_weight(),
        SymEngine::add(
            SymEngine::mul({
                weight,
                make_exc_density(D, Omega, 1),
                SymEngine::trace(
                    SymEngine::matrix_add({
                        SymEngine::matrix_mul({Omega_ab, D}),
                        SymEngine::matrix_mul({Omega_a, D_b})
                    })
                )
            }),
            SymEngine::mul({
                weight,
                make_exc_density(D, Omega, 2),
                SymEngine::trace(SymEngine::matrix_mul({Omega_a, D})),
                SymEngine::trace(SymEngine::matrix_add({
                    SymEngine::matrix_mul({Omega_b, D}),
                    SymEngine::matrix_mul({Omega, D_b})
                }))
            })
        )
    );
    auto hnuc_ab = SymEngine::rcp_dynamic_cast<const NonElecFunction>(hnuc_a->diff(b));
    // Equation (207) without T and S^{a}W matrices, J. Chem. Phys. 129, 214108 (2008)
    auto response = remove_if(E_ab, SymEngine::set_basic({D_a, D_ab}));
    REQUIRE(SymEngine::eq(
        *response,
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({h_a, D_b})),
            SymEngine::trace(SymEngine::matrix_mul({h_ab, D})),
            SymEngine::trace(SymEngine::matrix_mul({V_a, D_b})),
            SymEngine::trace(SymEngine::matrix_mul({V_ab, D})),
            // `TwoElecEnergy` considers trace(G(N)*M) = trace(G(M)*N),
            // Equation (210), J. Chem. Phys. 129, 214108 (2008)
            make_2el_energy(Ga, D_b),
            make_2el_energy(Gab, D),
            make_2el_energy(Ga_Db, D),
            Exc_ab_D,
            hnuc_ab
        })
    ));
}

TEST_CASE("Test ReplaceVisitor and replace()", "[ReplaceVisitor]")
{
    auto b = make_perturbation(std::string("b"));
    auto c = make_perturbation(std::string("c"));
    auto dependencies = PertDependency({std::make_pair(b, 99), std::make_pair(c, 99)});
    auto D = make_1el_density(std::string("D"));
    auto h = make_1el_operator(std::string("h"), dependencies);
    auto V = make_1el_operator(std::string("V"), dependencies);
    auto G = make_2el_operator(std::string("G"), D, dependencies);
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(std::string("Omega"), dependencies);
    auto Fxc = make_xc_potential(std::string("Fxc"), D, Omega, weight);
    auto F_ks = SymEngine::matrix_add({h, G, V, Fxc});
    auto T = make_t_matrix(dependencies);
    // Equation (94), J. Chem. Phys. 129, 214108 (2008)
    auto F = SymEngine::matrix_add({h, G, V, Fxc, T});
    auto S = make_1el_operator(std::string("S"), dependencies);
    // Equation (229), J. Chem. Phys. 129, 214108 (2008)
    auto Y = make_tdscf_equation(F, D, S);
    auto Dt = make_dt_operator(D);
    auto St = make_dt_operator(S);
    auto Y_b = Y->diff(b);
    auto DP = make_1el_density(std::string("D_{P}"));
    auto D_b = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(b));
    auto DP_b = SymEngine::rcp_dynamic_cast<const ElectronicState>(DP->diff(b));
    // Equation (165) without the frequency, J. Chem. Phys. 129, 214108 (2008)
    auto M_b = replace(
        // Unperturbed `Dt`, `St` and `T` will give zero frequency
        remove_if(Y_b, SymEngine::set_basic({Dt, St, T})),
        SymEngine::map_basic_basic({{D_b, DP_b}})
    );
    auto h_b = h->diff(b);
    auto V_b = V->diff(b);
    auto Gb = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({b})
    );
    auto G_DPb = make_2el_operator(G->get_name(), DP_b, dependencies);
    auto Omega_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(b));
    auto Fxc_b = Fxc->diff(b);
    // Equation (A26), J. Chem. Phys. 129, 214108 (2008)
    auto Fxc_DPb = replace(Fxc_b, SymEngine::map_basic_basic({{D_b, DP_b}}));
    REQUIRE(SymEngine::eq(
        *Fxc_DPb,
        // Equation (57) by replacing the first order perturbed density matrix
        // with the particular part, J. Chem. Phys. 140, 034103 (2014)
        *SymEngine::make_rcp<const ExchCorrPotential>(
            Fxc->get_name(),
            Fxc->get_state(),
            Fxc->get_overlap_distribution(),
            Fxc->get_weight(),
            SymEngine::matrix_add({
                SymEngine::matrix_mul({
                    SymEngine::mul({
                        weight,
                        make_exc_density(D, Omega, 2),
                        SymEngine::add(
                            make_density_vector(D, Omega_b),
                            make_density_vector(DP_b, Omega)
                        )
                    }),
                    Omega
                }),
                SymEngine::matrix_mul({
                    SymEngine::mul(weight, make_exc_density(D, Omega, 1)), Omega_b
                })
            })
        )
    ));
    auto T_b = T->diff(b);
    auto S_b = S->diff(b);
    auto DPt = make_dt_operator(DP);
    auto DPt_b = DPt->diff(b);
    auto St_b = St->diff(b);
    auto Fb = SymEngine::matrix_add({h_b, V_b, Gb, G_DPb, Fxc_DPb, T_b});
    auto minus_one_half = SymEngine::div(SymEngine::minus_one, SymEngine::two);
    REQUIRE(SymEngine::eq(
        *M_b,
        *SymEngine::matrix_add({
            SymEngine::matrix_mul({Fb, D, S}),
            SymEngine::matrix_mul({F_ks, DP_b, S}),
            SymEngine::matrix_mul({F_ks, D, S_b}),
            SymEngine::matrix_mul({SymEngine::minus_one, S, D, Fb}),
            SymEngine::matrix_mul({SymEngine::minus_one, S, DP_b, F_ks}),
            SymEngine::matrix_mul({SymEngine::minus_one, S_b, D, F_ks}),
            SymEngine::matrix_mul({minus_one_half, S, DPt_b, S}),
            SymEngine::matrix_mul({minus_one_half, S, D, St_b}),
            SymEngine::matrix_mul({minus_one_half, S, DPt_b, S}),
            SymEngine::matrix_mul({minus_one_half, St_b, D, S})
        })
    ));
    // A simple test for replacing more than one symbols
    auto GP = make_2el_operator(G->get_name(), DP, dependencies);
    auto FP_xc = make_xc_potential(Fxc->get_name(), DP, Omega, weight);
    auto FP_ks = SymEngine::matrix_add({h, GP, V, FP_xc});
    auto FP = SymEngine::matrix_add({h, GP, V, FP_xc, T});
    auto YP = make_tdscf_equation(FP, DP, S);
    REQUIRE(SymEngine::eq(
        *replace(
            remove_if(Y_b, SymEngine::set_basic({Dt, St, T})),
            SymEngine::map_basic_basic({{D, DP}, {D_b, DP_b}})
        ),
        *remove_if(YP->diff(b), SymEngine::set_basic({DPt, St, T}))
    ));
    auto Y_bc = Y_b->diff(c);
    auto D_c = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(c));
    auto D_bc = D_b->diff(c);
    auto DP_bc = SymEngine::rcp_dynamic_cast<const ElectronicState>(DP_b->diff(c));
    // Equation (187) without the frequency, J. Chem. Phys. 129, 214108 (2008)
    auto M_bc = replace(
        remove_if(Y_bc, SymEngine::set_basic({Dt, St, T})),
        SymEngine::map_basic_basic({{D_bc, DP_bc}})
    );
    auto h_bc = h_b->diff(c);
    auto V_bc = V_b->diff(c);
    auto Gbc = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({b, c})
    );
    auto Gb_Dc = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_c, dependencies, SymEngine::multiset_basic({b})
    );
    auto Gc_Db = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_b, dependencies, SymEngine::multiset_basic({c})
    );
    auto G_DPbc = make_2el_operator(G->get_name(), DP_bc, dependencies);
    auto Fxc_bc = Fxc_b->diff(c);
    // The second term of Equation (A31), J. Chem. Phys. 129, 214108 (2008)
    auto Fxc_DPbc = replace(Fxc_bc, SymEngine::map_basic_basic({{D_bc, DP_bc}}));
    auto T_bc = T_b->diff(c);
    auto S_c = S->diff(c);
    auto S_bc = S_b->diff(c);
    auto Dt_b = Dt->diff(b);
    auto Dt_c = Dt->diff(c);
    auto DPt_bc = DPt_b->diff(c);
    auto St_c = St->diff(c);
    auto St_bc = St_b->diff(c);
    auto Fbc = SymEngine::matrix_add({
        h_bc, V_bc, Gbc, Gb_Dc, Gc_Db, G_DPbc, Fxc_DPbc, T_bc
    });
    auto F_b = F->diff(b);
    auto F_c = F->diff(c);
    REQUIRE(SymEngine::eq(
        *M_bc,
        *SymEngine::matrix_add({
            SymEngine::matrix_mul({Fbc, D, S}),
            SymEngine::matrix_mul({F_ks, DP_bc, S}),
            SymEngine::matrix_mul({F_b, D_c, S}),
            SymEngine::matrix_mul({F_b, D, S_c}),
            SymEngine::matrix_mul({F_c, D_b, S}),
            SymEngine::matrix_mul({F_c, D, S_b}),
            SymEngine::matrix_mul({F_ks, D_b, S_c}),
            SymEngine::matrix_mul({F_ks, D_c, S_b}),
            SymEngine::matrix_mul({F_ks, D, S_bc}),
            SymEngine::matrix_mul({SymEngine::minus_one, S, D, Fbc}),
            SymEngine::matrix_mul({SymEngine::minus_one, S, DP_bc, F_ks}),
            SymEngine::matrix_mul({SymEngine::minus_one, S_c, D, F_b}),
            SymEngine::matrix_mul({SymEngine::minus_one, S, D_c, F_b}),
            SymEngine::matrix_mul({SymEngine::minus_one, S_b, D, F_c}),
            SymEngine::matrix_mul({SymEngine::minus_one, S, D_b, F_c}),
            SymEngine::matrix_mul({SymEngine::minus_one, S_c, D_b, F_ks}),
            SymEngine::matrix_mul({SymEngine::minus_one, S_b, D_c, F_ks}),
            SymEngine::matrix_mul({SymEngine::minus_one, S_bc, D, F_ks}),
            SymEngine::matrix_mul({minus_one_half, S_b, Dt_c, S}),
            SymEngine::matrix_mul({minus_one_half, S_b, D, St_c}),
            SymEngine::matrix_mul({minus_one_half, S_c, Dt_b, S}),
            SymEngine::matrix_mul({minus_one_half, S_c, D, St_b}),
            SymEngine::matrix_mul({minus_one_half, S, DPt_bc, S}),
            SymEngine::matrix_mul({minus_one_half, S, Dt_b, S_c}),
            SymEngine::matrix_mul({minus_one_half, S, D_b, St_c}),
            SymEngine::matrix_mul({minus_one_half, S, D_c, St_b}),
            SymEngine::matrix_mul({minus_one_half, S, Dt_c, S_b}),
            SymEngine::matrix_mul({minus_one_half, S, D, St_bc}),
            SymEngine::matrix_mul({minus_one_half, S, Dt_c, S_b}),
            SymEngine::matrix_mul({minus_one_half, St_c, D, S_b}),
            SymEngine::matrix_mul({minus_one_half, S, Dt_b, S_c}),
            SymEngine::matrix_mul({minus_one_half, St_b, D, S_c}),
            SymEngine::matrix_mul({minus_one_half, S, DPt_bc, S}),
            SymEngine::matrix_mul({minus_one_half, S_c, Dt_b, S}),
            SymEngine::matrix_mul({minus_one_half, St_c, D_b, S}),
            SymEngine::matrix_mul({minus_one_half, St_b, D_c, S}),
            SymEngine::matrix_mul({minus_one_half, S_b, Dt_c, S}),
            SymEngine::matrix_mul({minus_one_half, St_bc, D, S})
        })
    ));
}

TEST_CASE("Test FindAllVisitor and find_all()", "[FindAllVisitor]")
{
    auto b = make_perturbation(std::string("b"));
    auto c = make_perturbation(std::string("c"));
    auto dependencies = PertDependency({std::make_pair(b, 99), std::make_pair(c, 99)});
    auto D = make_1el_density(std::string("D"));
    auto h = make_1el_operator(std::string("h"), dependencies);
    auto V = make_1el_operator(std::string("V"), dependencies);
    auto G = make_2el_operator(std::string("G"), D, dependencies);
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(std::string("Omega"), dependencies);
    auto Exc = make_xc_energy(std::string("Exc"), D, Omega, weight);
    auto Fxc = make_xc_potential(std::string("Fxc"), D, Omega, weight);
    auto T = make_t_matrix(dependencies);
    auto hnuc = make_nonel_function(std::string("hnuc"), dependencies);
    auto E = make_ks_energy(h, V, G, D, Exc, hnuc);
    // Equation (94), J. Chem. Phys. 129, 214108 (2008)
    auto F = SymEngine::matrix_add({h, G, V, Fxc, T});
    auto S = make_1el_operator(std::string("S"), dependencies);

    // Test finding `TwoElecOperator` and `TwoElecEnergy` objects
    REQUIRE(is_equal(
        find_all(E, G), FindAllResult({{0, SymEngine::set_basic({G})}})
    ));
    REQUIRE(is_equal(
        find_all(E, make_2el_energy(G)),
        FindAllResult({{0, SymEngine::set_basic({make_2el_energy(G)})}})
    ));
    auto Gb = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({b})
    );
    REQUIRE(is_equal(
        find_all(E->diff(b), G),
        FindAllResult({{0, SymEngine::set_basic({G})}, {1, SymEngine::set_basic({Gb})}})
    ));

    // Equation (229), J. Chem. Phys. 129, 214108 (2008)
    auto Y = make_tdscf_equation(F, D, S);
    // (1) The first order
    auto Y_b = Y->diff(b);
    auto D_b = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(b));
    REQUIRE(is_equal(
        find_all(Y_b, D),
        FindAllResult({{0, SymEngine::set_basic({D})}, {1, SymEngine::set_basic({D_b})}})
    ));
    auto h_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(h->diff(b));
    REQUIRE(is_equal(
        find_all(Y_b, h),
        FindAllResult({{0, SymEngine::set_basic({h})}, {1, SymEngine::set_basic({h_b})}})
    ));
    auto V_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(V->diff(b));
    REQUIRE(is_equal(
        find_all(Y_b, V),
        FindAllResult({{0, SymEngine::set_basic({V})}, {1, SymEngine::set_basic({V_b})}})
    ));
    auto G_Db = make_2el_operator(G->get_name(), D_b, dependencies);
    REQUIRE(is_equal(
        find_all(Y_b, G),
        FindAllResult({
            {0, SymEngine::set_basic({G, G_Db})}, {1, SymEngine::set_basic({Gb})}
        })
    ));
    REQUIRE(is_equal(
        find_all(Y_b, weight), FindAllResult({{0, SymEngine::set_basic({weight})}})
    ));
    auto Omega_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(b));
    REQUIRE(is_equal(
        find_all(Y_b, Omega),
        FindAllResult({
            {0, SymEngine::set_basic({Omega})}, {1, SymEngine::set_basic({Omega_b})}
        })
    ));
    auto Fxc_b = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Fxc->diff(b));
    REQUIRE(is_equal(
        find_all(Y_b, Fxc),
        FindAllResult({
            {0, SymEngine::set_basic({Fxc})}, {1, SymEngine::set_basic({Fxc_b})}
        })
    ));
    auto T_b = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(T->diff(b));
    REQUIRE(is_equal(
        find_all(Y_b, T),
        FindAllResult({{0, SymEngine::set_basic({T})}, {1, SymEngine::set_basic({T_b})}})
    ));
    auto S_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(S->diff(b));
    REQUIRE(is_equal(
        find_all(Y_b, S),
        FindAllResult({{0, SymEngine::set_basic({S})}, {1, SymEngine::set_basic({S_b})}})
    ));
    // (2) The second order
    auto Y_bc = Y_b->diff(c);
    auto D_c = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(c));
    auto D_bc = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_b->diff(c));
    REQUIRE(is_equal(
        find_all(Y_bc, D),
        FindAllResult({
            {0, SymEngine::set_basic({D})},
            {1, SymEngine::set_basic({D_b, D_c})},
            {2, SymEngine::set_basic({D_bc})}
        })
    ));
    auto h_c = SymEngine::rcp_dynamic_cast<const OneElecOperator>(h->diff(c));
    auto h_bc = SymEngine::rcp_dynamic_cast<const OneElecOperator>(h_b->diff(c));
    REQUIRE(is_equal(
        find_all(Y_bc, h),
        FindAllResult({
            {0, SymEngine::set_basic({h})},
            {1, SymEngine::set_basic({h_b, h_c})},
            {2, SymEngine::set_basic({h_bc})}
        })
    ));
    auto V_c = SymEngine::rcp_dynamic_cast<const OneElecOperator>(V->diff(c));
    auto V_bc = SymEngine::rcp_dynamic_cast<const OneElecOperator>(V_b->diff(c));
    REQUIRE(is_equal(
        find_all(Y_bc, V),
        FindAllResult({
            {0, SymEngine::set_basic({V})},
            {1, SymEngine::set_basic({V_b, V_c})},
            {2, SymEngine::set_basic({V_bc})}
        })
    ));
    auto Gc = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({c})
    );
    auto G_Dc = make_2el_operator(G->get_name(), D_c, dependencies);
    auto Gb_Dc = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_c, dependencies, SymEngine::multiset_basic({b})
    );
    auto Gc_Db = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D_b, dependencies, SymEngine::multiset_basic({c})
    );
    auto Gbc = SymEngine::make_rcp<const TwoElecOperator>(
        G->get_name(), D, dependencies, SymEngine::multiset_basic({b, c})
    );
    auto G_Dbc = make_2el_operator(G->get_name(), D_bc, dependencies);
    REQUIRE(is_equal(
        find_all(Y_bc, G),
        FindAllResult({
            {0, SymEngine::set_basic({G, G_Db, G_Dc, G_Dbc})},
            {1, SymEngine::set_basic({Gb, Gc, Gb_Dc, Gc_Db})},
            {2, SymEngine::set_basic({Gbc})}
        })
    ));
    REQUIRE(is_equal(
        find_all(Y_bc, weight), FindAllResult({{0, SymEngine::set_basic({weight})}})
    ));
    auto Omega_c = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(c));
    auto Omega_bc = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega_b->diff(c));
    REQUIRE(is_equal(
        find_all(Y_bc, Omega),
        FindAllResult({
            {0, SymEngine::set_basic({Omega})},
            {1, SymEngine::set_basic({Omega_b, Omega_c})},
            {2, SymEngine::set_basic({Omega_bc})}
        })
    ));
    auto Fxc_c = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Fxc->diff(c));
    auto Fxc_bc = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Fxc_b->diff(c));
    REQUIRE(is_equal(
        find_all(Y_bc, Fxc),
        FindAllResult({
            {0, SymEngine::set_basic({Fxc})},
            {1, SymEngine::set_basic({Fxc_b, Fxc_c})},
            {2, SymEngine::set_basic({Fxc_bc})}
        })
    ));
    auto T_c = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(T->diff(c));
    auto T_bc = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(T_b->diff(c));
    REQUIRE(is_equal(
        find_all(Y_bc, T),
        FindAllResult({
            {0, SymEngine::set_basic({T})},
            {1, SymEngine::set_basic({T_b, T_c})},
            {2, SymEngine::set_basic({T_bc})}
        })
    ));
    auto S_c = SymEngine::rcp_dynamic_cast<const OneElecOperator>(S->diff(c));
    auto S_bc = SymEngine::rcp_dynamic_cast<const OneElecOperator>(S_b->diff(c));
    REQUIRE(is_equal(
        find_all(Y_bc, S),
        FindAllResult({
            {0, SymEngine::set_basic({S})},
            {1, SymEngine::set_basic({S_b, S_c})},
            {2, SymEngine::set_basic({S_bc})}
        })
    ));
}

//TEST_CASE("Test StringifyVisitor and stringify()", "[StringifyVisitor]")
//{
//}
