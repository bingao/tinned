#define CATCH_CONFIG_MAIN

#include <iostream>
#include <string>

#include <catch2/catch.hpp>

#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>

#include "Tinned.hpp"

using namespace Tinned;

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
    auto E = SymEngine::add(SymEngine::vec_basic({
        SymEngine::trace(SymEngine::matrix_mul(SymEngine::vec_basic({h, D}))),
        SymEngine::trace(SymEngine::matrix_mul(SymEngine::vec_basic({V, D}))),
        SymEngine::div(
            SymEngine::trace(SymEngine::matrix_mul(SymEngine::vec_basic({G, D}))),
            SymEngine::two
        ),
        Exc,
        hnuc
    }));
    std::cout << "E = " << stringify(E) << "\n\n";
    // Equation (81), J. Chem. Phys. 129, 214108 (2008)
    auto E_a = E->diff(a);
    std::cout << "E^{a} = " << stringify(E_a) << "\n\n";
    auto D_a = D->diff(a);
    auto E_0a = remove_if(E_a, SymEngine::set_basic({D_a}));
    std::cout << "E^{0,a} = " << stringify(E_0a) << "\n\n";
}

TEST_CASE("Test ReplaceVisitor and replace()", "[ReplaceVisitor]")
{

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
    auto Fxc = make_xc_potential(std::string("GGA"), D, Omega, weight);
    auto T = make_t_matrix(dependencies);
    // Equation (94), J. Chem. Phys. 129, 214108 (2008)
    auto F = SymEngine::matrix_add(SymEngine::vec_basic({h, G, V, Fxc, T}));
    std::cout << "F = " << stringify(F) << "\n\n";
    auto S = make_1el_operator(std::string("S"), dependencies);
    // Equation (229), J. Chem. Phys. 129, 214108 (2008)
    //auto Y = ;
}

//TEST_CASE("Test StringifyVisitor and stringify()", "[StringifyVisitor]")
//{
//}
