#define CATCH_CONFIG_MAIN

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <iostream>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/number.h>
#include <symengine/integer.h>
#include <symengine/constants.h>
#include <symengine/derivative.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/matrices/immutable_dense_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"
#include "TwoLevelAtom.hpp"

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

// unperturbed Hamiltonian [[E1, 0], [0, E2]]
inline SymEngine::RCP<const SymEngine::MatrixExpr> make_unperturbed_hamiltonian(
    const SymEngine::RCP<const SymEngine::Basic>& E1,
    const SymEngine::RCP<const SymEngine::Basic>& E2
)
{
    return SymEngine::immutable_dense_matrix(
        2, 2, {E1, SymEngine::zero, SymEngine::zero, E2}
    );
}

TEST_CASE("Test two-level system", "[FunctionEvaluator] and [OperatorEvaluator]")
{
    auto a = make_perturbation(std::string("a"));
    auto b = make_perturbation(std::string("b"));
    auto c = make_perturbation(std::string("c"));
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

    auto evaluator = TwoLevelFunction(SymEngine::vec_basic({a, b, c}));
    std::cout << stringify(evaluator.apply(E)) << "\n";
}
