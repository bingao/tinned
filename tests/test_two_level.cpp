#define CATCH_CONFIG_MAIN

#include <map>
#include <string>
#include <utility>

#include <iostream>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/real_double.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/symbol.h>
#include <symengine/matrices/immutable_dense_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"
#include "Tinned/TwoLevelAtom.hpp"

using namespace Tinned;

// This test adopts a two-state model system, from section 5.10.3, "Principles
// and Practices of Molecular Properties: Theory, Modeling and Simulations",
// Patrick Norman, Kenneth Ruud and Trond Saue.

// Ground state [[1, 0], [0, 0]]
inline SymEngine::RCP<const SymEngine::MatrixExpr> make_ground_state_density()
{
    return SymEngine::immutable_dense_matrix(
        2, 2, {SymEngine::one, SymEngine::zero, SymEngine::zero, SymEngine::zero}
    );
}

// Unperturbed Hamiltonian [[E1, 0], [0, E2]]
inline SymEngine::RCP<const SymEngine::MatrixExpr> make_unperturbed_hamiltonian(
    const SymEngine::RCP<const SymEngine::Basic>& E1,
    const SymEngine::RCP<const SymEngine::Basic>& E2
)
{
    return SymEngine::diagonal_matrix({E1, E2});
}

// Field operator [[0, mu], [mu, 0]], and mu must be real because we require
// that all field operators commute
inline SymEngine::RCP<const SymEngine::MatrixExpr> make_field_operator(
    const SymEngine::RCP<const SymEngine::Basic>& mu
)
{
    return SymEngine::immutable_dense_matrix(
        2, 2, {SymEngine::zero, mu, mu, SymEngine::zero}
    );
}

TEST_CASE("Test two-level system", "[FunctionEvaluator] and [OperatorEvaluator]")
{
    // Zero frequency is not allowed, see Equation (66) of ...
    //auto omega_a = SymEngine::real_double(0.1);
    //auto omega_b = SymEngine::real_double(0.2);
    //auto omega_c = SymEngine::real_double(0.3);
    auto omega_a = SymEngine::symbol("omega_a");
    auto omega_b = SymEngine::symbol("omega_b");
    auto omega_c = SymEngine::symbol("omega_c");
    auto a = make_perturbation(std::string("a"), omega_a);
    auto b = make_perturbation(std::string("b"), omega_b);
    auto c = make_perturbation(std::string("c"), omega_c);
    auto D = make_1el_density(std::string("D"));
    auto H0 = make_1el_operator(std::string("H0"));
    auto Va = make_1el_operator(
        std::string("Va"), PertDependency({std::make_pair(a, 1)})
    );
    auto Vb = make_1el_operator(
        std::string("Vb"), PertDependency({std::make_pair(b, 1)})
    );
    auto Vc = make_1el_operator(
        std::string("Vc"), PertDependency({std::make_pair(c, 1)})
    );
    auto E = SymEngine::add({
        SymEngine::trace(SymEngine::matrix_mul({H0, D})),
        SymEngine::trace(SymEngine::matrix_mul({Va, D})),
        SymEngine::trace(SymEngine::matrix_mul({Vb, D})),
        SymEngine::trace(SymEngine::matrix_mul({Vc, D}))
    });

    auto E1 = SymEngine::symbol("E_1");
    auto E2 = SymEngine::symbol("E_2");
    auto val_H0 = make_unperturbed_hamiltonian(E1, E2);
    auto mu_a = SymEngine::symbol("mu_a");
    auto mu_b = SymEngine::symbol("mu_b");
    auto mu_c = SymEngine::symbol("mu_c");
    auto val_Va = make_field_operator(mu_a);
    auto val_Vb = make_field_operator(mu_b);
    auto val_Vc = make_field_operator(mu_c);

    auto oper_eval = TwoLevelOperator(
        std::make_pair(H0, val_H0),
        std::map<SymEngine::RCP<const OneElecOperator>,
                 SymEngine::RCP<const SymEngine::MatrixExpr>,
                 SymEngine::RCPBasicKeyLess>({{Va, val_Va}, {Vb, val_Vb}, {Vc, val_Vc}}),
        std::make_pair(D, make_ground_state_density())
    );

    REQUIRE(SymEngine::eq(*oper_eval.apply(H0), *val_H0));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Va), *SymEngine::matrix_mul({a, val_Va})));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Vb), *SymEngine::matrix_mul({b, val_Vb})));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Vc), *SymEngine::matrix_mul({c, val_Vc})));
    REQUIRE(SymEngine::eq(*oper_eval.apply(D), *make_ground_state_density()));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Va->diff(a)), *val_Va));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Vb->diff(b)), *val_Vb));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Vc->diff(c)), *val_Vc));
    // D^{a} = [[0, mu_a/(E_1-E_2-omega_a)], [mu_a/(E_1-E_2+omega_a), 0]]
    REQUIRE(SymEngine::eq(
        *oper_eval.apply(D->diff(a)),
        *SymEngine::immutable_dense_matrix(
            2, 2,
            {
                SymEngine::zero,
                SymEngine::div(
                    SymEngine::mul(SymEngine::minus_one, mu_a),
                    SymEngine::sub(omega_a, SymEngine::sub(E1, E2))
                ),
                SymEngine::div(mu_a, SymEngine::sub(omega_a, SymEngine::sub(E2, E1))),
                SymEngine::zero
            }
        )
    ));
    // D^{ab} = [[mu_a*mu_b/(omega_a+omega_b)*(1/(E_1-E2+omega_a)+1/(E_1-E_2+omega_b)), 0],
    //           [0, mu_a*mu_b/(omega_a+omega_b)*(1/(E_1-E2-omega_a)+1/(E_1-E_2-omega_b))]]
//    REQUIRE(SymEngine::eq(
//        *oper_eval.apply((D->diff(a))->diff(b)),
//
//    ));

std::cout << "D^{ab} = " << stringify(oper_eval.apply((D->diff(a))->diff(b))) << "\n\n";

    auto fun_eval = TwoLevelFunction(
        std::make_pair(H0, val_H0),
        std::map<SymEngine::RCP<const OneElecOperator>,
                 SymEngine::RCP<const SymEngine::MatrixExpr>,
                 SymEngine::RCPBasicKeyLess>({{Va, val_Va}, {Vb, val_Vb}, {Vc, val_Vc}}),
        std::make_pair(D, make_ground_state_density())
    );

    std::cout << "E = " << stringify(E) << "\n";
    std::cout << "E = " << stringify(fun_eval.apply(E)) << "\n";
    std::cout << "E^a = " << stringify(E->diff(a)) << "\n";
    std::cout << "E^a = " << stringify(fun_eval.apply(E->diff(a))) << "\n";
}
