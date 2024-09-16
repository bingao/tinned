#define CATCH_CONFIG_MAIN

#include <map>
#include <string>
#include <utility>

#include <iostream>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/constants.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/symbol.h>
#include <symengine/functions.h>
#include <symengine/matrices/diagonal_matrix.h>
#include <symengine/matrices/immutable_dense_matrix.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/hadamard_product.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"
#include "Tinned/TwoLevelAtom.hpp"

using namespace Tinned;

// This test adopts a two-state model system, from section 5.10.3, "Principles
// and Practices of Molecular Properties: Theory, Modeling and Simulations",
// Patrick Norman, Kenneth Ruud and Trond Saue.

// Density matrix for a pure state
// [[rho, exp(-i*phi)*sqrt(rho*(1-rho))], [exp(i*phi)*sqrt(rho*(1-rho)), 1-rho]]
// See ...
inline SymEngine::RCP<const SymEngine::MatrixExpr> pure_state_density(
    const SymEngine::RCP<const SymEngine::Basic>& rho,
    const SymEngine::RCP<const SymEngine::Basic>& phi
)
{
    auto rho_22 = SymEngine::sub(SymEngine::one, rho);
    if (SymEngine::eq(*rho, *SymEngine::zero) ||
        SymEngine::eq(*rho_22, *SymEngine::zero)) {
        return SymEngine::diagonal_matrix({rho, rho_22});
    }
    else {
        auto rho_12 = SymEngine::sqrt(SymEngine::mul(rho, rho_22));
        if (SymEngine::eq(*phi, *SymEngine::zero)) {
            return SymEngine::immutable_dense_matrix(
                2, 2, {rho, rho_12, rho_12, rho_22}
            );
        }
        else {
            return SymEngine::immutable_dense_matrix(
                2, 2,
                {
                    rho,
                    SymEngine::mul(
                        SymEngine::sub(
                            SymEngine::cos(phi),
                            SymEngine::mul(SymEngine::I, SymEngine::sin(phi))
                        ),
                        rho_12
                    ),
                    SymEngine::mul(
                        SymEngine::add(
                            SymEngine::cos(phi),
                            SymEngine::mul(SymEngine::I, SymEngine::sin(phi))
                        ),
                        rho_12
                    ),
                    rho_22
                }
            );
        }
    }
}

// Density matrix for mixed states
// [[rho_11, rho_12], [conj(rho_12), 1-rho_11]]
// See ...
//inline SymEngine::RCP<const SymEngine::MatrixExpr> mixed_states_density(
//    const SymEngine::RCP<const SymEngine::Basic>& rho_11,
//    const SymEngine::RCP<const SymEngine::Basic>& rho_12
//)
//{
//    auto rho_22 = SymEngine::sub(SymEngine::one, rho_11);
//    return SymEngine::immutable_dense_matrix(
//        2, 2, {rho_11, rho_12, SymEngine::conjugate(rho_12), rho_22}
//    );
//}

// Unperturbed Hamiltonian [[E1, 0], [0, E2]]
inline SymEngine::RCP<const SymEngine::MatrixExpr> make_unperturbed_hamiltonian(
    const SymEngine::RCP<const SymEngine::Basic>& E1,
    const SymEngine::RCP<const SymEngine::Basic>& E2
)
{
    return SymEngine::diagonal_matrix({E1, E2});
}

// Field operator [[0, mu], [mu, 0]], and mu must be real because we require
// that any pair of field operators commute
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
    auto val_Ba = make_field_operator(mu_a);
    auto val_Bb = make_field_operator(mu_b);
    auto val_Bc = make_field_operator(mu_c);

    // Test ground state [[1, 0], [0, 0]]
    auto val_D = pure_state_density(SymEngine::one, SymEngine::zero);

    auto oper_eval = TwoLevelOperator(
        std::make_pair(H0, val_H0),
        std::map<SymEngine::RCP<const OneElecOperator>,
                 SymEngine::RCP<const SymEngine::MatrixExpr>,
                 SymEngine::RCPBasicKeyLess>({{Va, val_Ba}, {Vb, val_Bb}, {Vc, val_Bc}}),
        std::make_pair(D, val_D)
    );

    REQUIRE(SymEngine::eq(*oper_eval.apply(H0), *val_H0));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Va), *SymEngine::matrix_mul({a, val_Ba})));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Vb), *SymEngine::matrix_mul({b, val_Bb})));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Vc), *SymEngine::matrix_mul({c, val_Bc})));
    REQUIRE(SymEngine::eq(*oper_eval.apply(D), *val_D));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Va->diff(a)), *val_Ba));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Vb->diff(b)), *val_Bb));
    REQUIRE(SymEngine::eq(*oper_eval.apply(Vc->diff(c)), *val_Bc));
    // D^{a} = [[0, -mu_a/(omega_a+(E_2-E_1))], [mu_a/(omega_a-(E_2-E_1)), 0]]
    auto E_12 = SymEngine::sub(E1, E2);
    auto E_21 = SymEngine::sub(E2, E1);
    REQUIRE(SymEngine::eq(
        *oper_eval.apply(D->diff(a)),
        *SymEngine::immutable_dense_matrix(
            2, 2,
            {
                SymEngine::zero,
                SymEngine::div(
                    SymEngine::mul(SymEngine::minus_one, mu_a),
                    SymEngine::sub(omega_a, E_12)
                ),
                SymEngine::div(mu_a, SymEngine::sub(omega_a, E_21)),
                SymEngine::zero
            }
        )
    ));
    // D^{ab} = mu_a*mu_b/(omega_a+omega_b)*[
    //     [1/(omega_a+(E_2-E_1))+1/(omega_a-(E_2-E_1))+1/(omega_b+(E_2-E_1))+1/(omega_b-(E_2-E_1)), 0],
    //     [0, -1/(omega_a+(E_2-E_1))-1/(omega_a-(E_2-E_1))-1/(omega_b+(E_2-E_1))-1/(omega_b-(E_2-E_1))]
    // ]
    REQUIRE(SymEngine::eq(
        *oper_eval.apply(D->diff(a)->diff(b)),
        *SymEngine::diagonal_matrix({
            SymEngine::add(
                SymEngine::div(
                    SymEngine::add(
                        SymEngine::div(
                            SymEngine::mul(mu_a, mu_b), SymEngine::sub(omega_a, E_12)
                        ),
                        SymEngine::div(
                            SymEngine::mul(mu_a, mu_b), SymEngine::sub(omega_a, E_21)
                        )
                    ),
                    SymEngine::add(omega_a, omega_b)
                ),
                SymEngine::div(
                    SymEngine::add(
                        SymEngine::div(
                            SymEngine::mul(mu_a, mu_b), SymEngine::sub(omega_b, E_12)
                        ),
                        SymEngine::div(
                            SymEngine::mul(mu_a, mu_b), SymEngine::sub(omega_b, E_21)
                        )
                    ),
                    SymEngine::add(omega_a, omega_b)
                )
            ),
            SymEngine::add(
                SymEngine::div(
                    SymEngine::add(
                        SymEngine::div(
                            SymEngine::mul({SymEngine::minus_one, mu_a, mu_b}),
                            SymEngine::sub(omega_a, E_12)
                        ),
                        SymEngine::div(
                            SymEngine::mul({SymEngine::minus_one, mu_a, mu_b}),
                            SymEngine::sub(omega_a, E_21)
                        )
                    ),
                    SymEngine::add(omega_a, omega_b)
                ),
                SymEngine::div(
                    SymEngine::add(
                        SymEngine::div(
                            SymEngine::mul({SymEngine::minus_one, mu_a, mu_b}),
                            SymEngine::sub(omega_b, E_12)
                        ),
                        SymEngine::div(
                            SymEngine::mul({SymEngine::minus_one, mu_a, mu_b}),
                            SymEngine::sub(omega_b, E_21)
                        )
                    ),
                    SymEngine::add(omega_a, omega_b)
                )
            )
        })
    ));
    // D^{ab} = D^{ba}
    REQUIRE(SymEngine::eq(
        *oper_eval.apply(D->diff(a)->diff(b)), *oper_eval.apply(D->diff(b)->diff(a))
    ));
    // D^{abc} = 2*[Vc^{c}, D^{ab}]\circ[(omega_a+omega_b+omega_c-(E_m-E_n))^{-1}]
    //         + 2*[Vb^{b}, D^{ac}]\circ[(omega_a+omega_b+omega_c-(E_m-E_n))^{-1}]
    //         + 2*[Va^{a}, D^{bc}]\circ[(omega_a+omega_b+omega_c-(E_m-E_n))^{-1}]
    auto val_Dab = oper_eval.apply(D->diff(a)->diff(b));
    auto val_Dac = oper_eval.apply(D->diff(a)->diff(c));
    auto val_Dbc = oper_eval.apply(D->diff(b)->diff(c));
    auto minus_identity = SymEngine::diagonal_matrix({
        SymEngine::minus_one, SymEngine::minus_one
    });
    auto val_fabc = SymEngine::immutable_dense_matrix(
        2, 2,
        {
            SymEngine::div(SymEngine::two, SymEngine::add({omega_a, omega_b, omega_c})),
            SymEngine::div(
                SymEngine::two,
                SymEngine::sub(SymEngine::add({omega_a, omega_b, omega_c}), E_12)
            ),
            SymEngine::div(
                SymEngine::two,
                SymEngine::sub(SymEngine::add({omega_a, omega_b, omega_c}), E_21)
            ),
            SymEngine::div(SymEngine::two, SymEngine::add({omega_a, omega_b, omega_c}))
        }
    );
    REQUIRE(SymEngine::eq(
        *oper_eval.apply(D->diff(a)->diff(b)->diff(c)),
        *SymEngine::matrix_add({
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Bc, val_Dab}),
                    // `minus_identity` makes sure `SymEngine::matrix_mul`
                    // returns only a 2x2 matrix
                    SymEngine::matrix_mul({minus_identity, val_Dab, val_Bc})
                }),
                val_fabc
            }),
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Bb, val_Dac}),
                    SymEngine::matrix_mul({minus_identity, val_Dac, val_Bb})
                }),
                val_fabc
            }),
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Ba, val_Dbc}),
                    SymEngine::matrix_mul({minus_identity, val_Dbc, val_Ba})
                }),
                val_fabc
            })
        })
    ));
    // D^{aab} = 4*[Va^{a}, D^{ab}]\circ[(2*omega_a+omega_b-(E_m-E_n))^{-1}]
    //         + 2*[Vb^{b}, D^{aa}]\circ[(2*omega_a+omega_b-(E_m-E_n))^{-1}]
    auto val_Daa = oper_eval.apply(D->diff(a)->diff(a));
    auto val_faba = SymEngine::immutable_dense_matrix(
        2, 2,
        {
            SymEngine::div(
                SymEngine::integer(4), SymEngine::add({omega_a, omega_a, omega_b})
            ),
            SymEngine::div(
                SymEngine::integer(4),
                SymEngine::sub(SymEngine::add({omega_a, omega_a, omega_b}), E_12)
            ),
            SymEngine::div(
                SymEngine::integer(4),
                SymEngine::sub(SymEngine::add({omega_a, omega_a, omega_b}), E_21)
            ),
            SymEngine::div(
                SymEngine::integer(4), SymEngine::add({omega_a, omega_a, omega_b})
            )
        }
    );
    auto val_faab = SymEngine::immutable_dense_matrix(
        2, 2,
        {
            SymEngine::div(SymEngine::two, SymEngine::add({omega_a, omega_a, omega_b})),
            SymEngine::div(
                SymEngine::two,
                SymEngine::sub(SymEngine::add({omega_a, omega_a, omega_b}), E_12)
            ),
            SymEngine::div(
                SymEngine::two,
                SymEngine::sub(SymEngine::add({omega_a, omega_a, omega_b}), E_21)
            ),
            SymEngine::div(SymEngine::two, SymEngine::add({omega_a, omega_a, omega_b}))
        }
    );
    REQUIRE(SymEngine::eq(
        *oper_eval.apply(D->diff(a)->diff(a)->diff(b)),
        *SymEngine::matrix_add({
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Ba, val_Dab}),
                    SymEngine::matrix_mul({minus_identity, val_Dab, val_Ba})
                }),
                val_faba
            }),
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Bb, val_Daa}),
                    SymEngine::matrix_mul({minus_identity, val_Daa, val_Bb})
                }),
                val_faab
            })
        })
    ));

    auto fun_eval = TwoLevelFunction(
        std::make_pair(H0, val_H0),
        std::map<SymEngine::RCP<const OneElecOperator>,
                 SymEngine::RCP<const SymEngine::MatrixExpr>,
                 SymEngine::RCPBasicKeyLess>({{Va, val_Ba}, {Vb, val_Bb}, {Vc, val_Bc}}),
        std::make_pair(D, val_D)
    );

    REQUIRE(SymEngine::eq(*fun_eval.apply(E), *E1));
    // E^{a} = tr(D^{a}*Va) + tr(D^{a}*Vb) + tr(D^{a}*Vc)
    REQUIRE(SymEngine::eq(
        *fun_eval.apply(E->diff(a)),
        *SymEngine::add({
            SymEngine::mul(
                a,
                SymEngine::add(
                    SymEngine::div(
                        SymEngine::mul(
                            SymEngine::minus_one, SymEngine::pow(mu_a, SymEngine::two)
                        ),
                        SymEngine::sub(omega_a, E_12)
                    ),
                    SymEngine::div(
                        SymEngine::pow(mu_a, SymEngine::two),
                        SymEngine::sub(omega_a, E_21)
                    )
                )
            ),
            SymEngine::mul(
                b,
                SymEngine::add(
                    SymEngine::div(
                        SymEngine::mul({SymEngine::minus_one, mu_a, mu_b}),
                        SymEngine::sub(omega_a, E_12)
                    ),
                    SymEngine::div(
                        SymEngine::mul(mu_a, mu_b), SymEngine::sub(omega_a, E_21)
                    )
                )
            ),
            SymEngine::mul(
                c,
                SymEngine::add(
                    SymEngine::div(
                        SymEngine::mul({SymEngine::minus_one, mu_a, mu_c}),
                        SymEngine::sub(omega_a, E_12)
                    ),
                    SymEngine::div(
                        SymEngine::mul(mu_a, mu_c), SymEngine::sub(omega_a, E_21)
                    )
                )
            )
        })
    ));
    // E^{ab} = tr(D^{b}*Va^{a}) + tr(D^{a}*Vb^{b}) + tr(D^{ab}*H0)
    auto val_Da = oper_eval.apply(D->diff(a));
    auto val_Db = oper_eval.apply(D->diff(b));
    REQUIRE(SymEngine::eq(
        *fun_eval.apply(E->diff(a)->diff(b)),
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({val_Db, val_Ba})),
            SymEngine::trace(SymEngine::matrix_mul({val_Da, val_Bb})),
            SymEngine::trace(SymEngine::matrix_mul({val_Dab, val_H0}))
        })
    ));
    // E^{abc} = tr(D^{abc}*Va) + tr(D^{abc}*Vb) + tr(D^{abc}*Vc)
    auto val_Dabc = oper_eval.apply(D->diff(a)->diff(b)->diff(c));
    REQUIRE(SymEngine::eq(
        *fun_eval.apply(E->diff(a)->diff(b)->diff(c)),
        *SymEngine::add({
            SymEngine::mul(
                a, SymEngine::trace(SymEngine::matrix_mul({val_Dabc, val_Ba}))
            ),
            SymEngine::mul(
                b, SymEngine::trace(SymEngine::matrix_mul({val_Dabc, val_Bb}))
            ),
            SymEngine::mul(
                c, SymEngine::trace(SymEngine::matrix_mul({val_Dabc, val_Bc}))
            )
        })
    ));

    //How to test mixed_states_density()? Will V's be different?
}
