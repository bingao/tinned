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
    return SymEngine::add(SymEngine::vec_basic({
        SymEngine::trace(SymEngine::matrix_mul(SymEngine::vec_basic({h, D}))),
        SymEngine::trace(SymEngine::matrix_mul(SymEngine::vec_basic({V, D}))),
        SymEngine::div(
            SymEngine::trace(SymEngine::matrix_mul(SymEngine::vec_basic({G, D}))),
            SymEngine::two
        ),
        Exc,
        hnuc
    }));
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
    return SymEngine::matrix_add(SymEngine::vec_basic({
        SymEngine::matrix_mul(SymEngine::vec_basic({F, D, S})),
        SymEngine::matrix_mul(SymEngine::vec_basic({SymEngine::minus_one, S, D, F})),
        SymEngine::matrix_mul(SymEngine::vec_basic({SymEngine::minus_one, S, Dt, S})),
        SymEngine::matrix_mul(SymEngine::vec_basic({minus_one_half, St, D, S})),
        SymEngine::matrix_mul(SymEngine::vec_basic({minus_one_half, S, D, St}))
    }));
}

TEST_CASE("Test KeepVisitor and keep_if()", "[KeepVisitor]")
{

}

TEST_CASE("Test RemoveVisitor and remove_if()", "[RemoveVisitor]")
{
    auto a = make_perturbation(std::string("a"));
    auto dependencies = PertDependency({std::make_pair(a, 99)});
    auto D = make_1el_density(std::string("D"));
    auto h = make_1el_operator(std::string("h"), dependencies);
    auto V = make_1el_operator(std::string("V"), dependencies);
    auto G = make_2el_operator(std::string("G"), D, dependencies);
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(std::string("Omega"), dependencies);
    auto Exc = make_xc_energy(std::string("GGA"), D, Omega, weight);
    auto hnuc = make_nonel_function(std::string("hnuc"), dependencies);
    // Equation (80), J. Chem. Phys. 129, 214108 (2008)
    auto E = make_ks_energy(h, V, G, D, Exc, hnuc);
    // Equation (81), J. Chem. Phys. 129, 214108 (2008)
    auto E_a = E->diff(a);
    auto D_a = D->diff(a);
    //auto E_0a = remove_if(E_a, SymEngine::set_basic({D_a}));
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
    auto Fxc = make_xc_potential(std::string("GGA"), D, Omega, weight);
    auto T = make_t_matrix(dependencies);
    // Equation (94), J. Chem. Phys. 129, 214108 (2008)
    auto F = SymEngine::matrix_add(SymEngine::vec_basic({h, G, V, Fxc, T}));
    auto S = make_1el_operator(std::string("S"), dependencies);
    // Equation (229), J. Chem. Phys. 129, 214108 (2008)
    auto Y = make_tdscf_equation(F, D, S);
    // Equation (165), J. Chem. Phys. 129, 214108 (2008)
    auto Y_b = Y->diff(b);
    auto DP = make_1el_density(std::string("D_{P}"));
    auto D_b = D->diff(b);
    auto DP_b = DP->diff(b);
    //auto M_b = replace(Y_b, SymEngine::map_basic_basic({{D_b, DP_b}}));
    // Equation (187), J. Chem. Phys. 129, 214108 (2008)
    auto Y_bc = Y_b->diff(c);
    auto D_bc = D_b->diff(c);
    auto DP_bc = DP_b->diff(c);
    //auto M_bc = replace(Y_bc, SymEngine::map_basic_basic({{D_bc, DP_bc}}));
}

TEST_CASE("Test FindAllVisitor and find_all()", "[FindAllVisitor]")
{
    auto b = make_perturbation(std::string("b"));
    auto c = make_perturbation(std::string("c"));
    auto dependencies = PertDependency({std::make_pair(b, 99), std::make_pair(c, 99)});
    auto D = make_1el_density(std::string("D"));
    auto h = make_1el_operator(std::string("h"), dependencies);
    auto V = make_1el_operator(std::string("V"), dependencies);
    auto G_name = std::string("G");
    auto G = make_2el_operator(G_name, D, dependencies);
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(std::string("Omega"), dependencies);
    auto Fxc = make_xc_potential(std::string("GGA"), D, Omega, weight);
    auto T = make_t_matrix(dependencies);
    // Equation (94), J. Chem. Phys. 129, 214108 (2008)
    auto F = SymEngine::matrix_add(SymEngine::vec_basic({h, G, V, Fxc, T}));
    auto S = make_1el_operator(std::string("S"), dependencies);
    // Equation (229), J. Chem. Phys. 129, 214108 (2008)
    auto Y = make_tdscf_equation(F, D, S);
    // (1) The first order
    auto Y_b = Y->diff(b);
    auto D_b = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(b));
    REQUIRE(SymEngine::unified_eq(
        find_all<const ElectronicState>(Y_b, D),
        SameTypeSet<const ElectronicState>({D, D_b})
    ));
    auto h_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(h->diff(b));
    REQUIRE(SymEngine::unified_eq(
        find_all<const OneElecOperator>(Y_b, h),
        SameTypeSet<const OneElecOperator>({h, h_b})
    ));
    auto V_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(V->diff(b));
    REQUIRE(SymEngine::unified_eq(
        find_all<const OneElecOperator>(Y_b, V),
        SameTypeSet<const OneElecOperator>({V, V_b})
    ));
    auto G_b = SymEngine::make_rcp<const TwoElecOperator>(
        G_name, D, dependencies, SymEngine::multiset_basic({b})
    );
    auto G_Db = make_2el_operator(G_name, D_b, dependencies);
    REQUIRE(SymEngine::unified_eq(
        find_all<const TwoElecOperator>(Y_b, G),
        SameTypeSet<const TwoElecOperator>({G, G_b, G_Db})
    ));
    REQUIRE(SymEngine::unified_eq(
        find_all<const NonElecFunction>(Y_b, weight),
        SameTypeSet<const NonElecFunction>({weight})
    ));
    auto Omega_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(Omega->diff(b));
    REQUIRE(SymEngine::unified_eq(
        find_all<const OneElecOperator>(Y_b, Omega),
        SameTypeSet<const OneElecOperator>({Omega, Omega_b})
    ));
    auto Fxc_b = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Fxc->diff(b));
    REQUIRE(SymEngine::unified_eq(
        find_all<const ExchCorrPotential>(Y_b, Fxc),
        SameTypeSet<const ExchCorrPotential>({Fxc, Fxc_b})
    ));
    auto T_b = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(T->diff(b));
    REQUIRE(SymEngine::unified_eq(
        find_all<const TemporumOverlap>(Y_b, T),
        SameTypeSet<const TemporumOverlap>({T, T_b})
    ));
    auto S_b = SymEngine::rcp_dynamic_cast<const OneElecOperator>(S->diff(b));
    REQUIRE(SymEngine::unified_eq(
        find_all<const OneElecOperator>(Y_b, S),
        SameTypeSet<const OneElecOperator>({S, S_b})
    ));
    // (2) The second order
    auto Y_bc = Y_b->diff(c);
    auto D_c = SymEngine::rcp_dynamic_cast<const ElectronicState>(D->diff(c));
    auto D_bc = SymEngine::rcp_dynamic_cast<const ElectronicState>(D_b->diff(c));
    REQUIRE(SymEngine::unified_eq(
        find_all<const ElectronicState>(Y_bc, D),
        SameTypeSet<const ElectronicState>({D, D_b, D_c, D_bc})
    ));
}

//TEST_CASE("Test StringifyVisitor and stringify()", "[StringifyVisitor]")
//{
//}
