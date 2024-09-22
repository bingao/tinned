#define CATCH_CONFIG_MAIN

#include <iostream>

#include <cstddef>
#include <map>
#include <string>
#include <utility>

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
    auto rho_11 = SymEngine::sub(SymEngine::one, rho);
    if (SymEngine::eq(*rho, *SymEngine::zero) ||
        SymEngine::eq(*rho_11, *SymEngine::zero)) {
        return SymEngine::diagonal_matrix({rho, rho_11});
    }
    else {
        auto rho_01 = SymEngine::sqrt(SymEngine::mul(rho, rho_11));
        if (SymEngine::eq(*phi, *SymEngine::zero)) {
            return SymEngine::immutable_dense_matrix(
                2, 2, {rho, rho_01, rho_01, rho_11}
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
                        rho_01
                    ),
                    SymEngine::mul(
                        SymEngine::add(
                            SymEngine::cos(phi),
                            SymEngine::mul(SymEngine::I, SymEngine::sin(phi))
                        ),
                        rho_01
                    ),
                    rho_11
                }
            );
        }
    }
}

// Density matrix for mixed states
// [[rho_11, rho_01], [conj(rho_01), 1-rho_11]]
// See ...
//inline SymEngine::RCP<const SymEngine::MatrixExpr> mixed_states_density(
//    const SymEngine::RCP<const SymEngine::Basic>& rho_11,
//    const SymEngine::RCP<const SymEngine::Basic>& rho_01
//)
//{
//    auto rho_11 = SymEngine::sub(SymEngine::one, rho_11);
//    return SymEngine::immutable_dense_matrix(
//        2, 2, {rho_11, rho_01, SymEngine::conjugate(rho_01), rho_11}
//    );
//}

// Unperturbed Hamiltonian [[E0, 0], [0, E1]]
inline SymEngine::RCP<const SymEngine::MatrixExpr> make_unperturbed_hamiltonian(
    const SymEngine::RCP<const SymEngine::Basic>& E0,
    const SymEngine::RCP<const SymEngine::Basic>& E1
)
{
    return SymEngine::diagonal_matrix({E0, E1});
}

// Field operator
inline SymEngine::RCP<const SymEngine::MatrixExpr> make_field_operator(
    const SymEngine::RCP<const SymEngine::Basic>& V_00,
    const SymEngine::RCP<const SymEngine::Basic>& V_01,
    const SymEngine::RCP<const SymEngine::Basic>& V_10,
    const SymEngine::RCP<const SymEngine::Basic>& V_11
)
{
    return SymEngine::immutable_dense_matrix(
        2, 2, {V_00, V_01, V_10, V_11}
    );
}

TEST_CASE("Test two-level system", "[TwoLevelFunction], [TwoLevelOperator]")
{
    // Zero frequency is not allowed, see Equation (66) of ...
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

    auto E0 = SymEngine::symbol("E_0");
    auto E1 = SymEngine::symbol("E_1");
    auto val_H0 = make_unperturbed_hamiltonian(E0, E1);
    auto Va_00 = SymEngine::symbol("Va_00");
    auto Va_01 = SymEngine::symbol("Va_01");
    auto Va_10 = SymEngine::symbol("Va_10");
    auto Va_11 = SymEngine::symbol("Va_11");
    auto Vb_00 = SymEngine::symbol("Vb_00");
    auto Vb_01 = SymEngine::symbol("Vb_01");
    auto Vb_10 = SymEngine::symbol("Vb_10");
    auto Vb_11 = SymEngine::symbol("Vb_11");
    auto Vc_00 = SymEngine::symbol("Vc_00");
    auto Vc_01 = SymEngine::symbol("Vc_01");
    auto Vc_10 = SymEngine::symbol("Vc_10");
    auto Vc_11 = SymEngine::symbol("Vc_11");
    auto val_Ba = make_field_operator(Va_00, Va_01, Va_10, Va_11);
    auto val_Bb = make_field_operator(Vb_00, Vb_01, Vb_10, Vb_11);
    auto val_Bc = make_field_operator(Vc_00, Vc_01, Vc_10, Vc_11);

    // We choose ground state [[1, 0], [0, 0]]
    auto val_D = pure_state_density(SymEngine::one, SymEngine::zero);

    // Test evaluator for operators
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
    // D^{a} = [[0, -Va_01/(omega_a+(E1-E0))], [Va_10/(omega_a-(E1-E0)), 0]]
    auto E_01 = SymEngine::sub(E0, E1);
    auto E_10 = SymEngine::sub(E1, E0);
    REQUIRE(SymEngine::eq(
        *oper_eval.apply(D->diff(a)),
        *SymEngine::immutable_dense_matrix(
            2, 2,
            {
                SymEngine::zero,
                SymEngine::div(
                    SymEngine::mul(SymEngine::minus_one, Va_01),
                    SymEngine::sub(omega_a, E_01)
                ),
                SymEngine::div(Va_10, SymEngine::sub(omega_a, E_10)),
                SymEngine::zero
            }
        )
    ));
    // D^{ab} = [Vb^{b}, D^{a}]\circ[(omega_a+omega_b-(E_m-E_n))^{-1}]
    //        + [Va^{a}, D^{b}]\circ[(omega_a+omega_b-(E_m-E_n))^{-1}]
    auto val_Da = oper_eval.apply(D->diff(a));
    auto val_Db = oper_eval.apply(D->diff(b));
    // Scalar matrix with the number as minus one
    auto minus_identity = SymEngine::diagonal_matrix({
        SymEngine::minus_one, SymEngine::minus_one
    });
    auto val_freq_denom = SymEngine::immutable_dense_matrix(
        2, 2,
        {
            SymEngine::div(SymEngine::one, SymEngine::add(omega_a, omega_b)),
            SymEngine::div(
                SymEngine::one,
                SymEngine::sub(SymEngine::add(omega_a, omega_b), E_01)
            ),
            SymEngine::div(
                SymEngine::one,
                SymEngine::sub(SymEngine::add(omega_a, omega_b), E_10)
            ),
            SymEngine::div(SymEngine::one, SymEngine::add(omega_a, omega_b))
        }
    );
    REQUIRE(SymEngine::eq(
        *oper_eval.apply(D->diff(a)->diff(b)),
        *SymEngine::matrix_add({
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Bb, val_Da}),
                    // `minus_identity` makes sure `SymEngine::matrix_mul`
                    // returns only a 2x2 matrix
                    SymEngine::matrix_mul({minus_identity, val_Da, val_Bb})
                }),
                val_freq_denom
            }),
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Ba, val_Db}),
                    SymEngine::matrix_mul({minus_identity, val_Db, val_Ba})
                }),
                val_freq_denom
            })
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
    val_freq_denom = SymEngine::immutable_dense_matrix(
        2, 2,
        {
            SymEngine::div(SymEngine::two, SymEngine::add({omega_a, omega_b, omega_c})),
            SymEngine::div(
                SymEngine::two,
                SymEngine::sub(SymEngine::add({omega_a, omega_b, omega_c}), E_01)
            ),
            SymEngine::div(
                SymEngine::two,
                SymEngine::sub(SymEngine::add({omega_a, omega_b, omega_c}), E_10)
            ),
            SymEngine::div(SymEngine::two, SymEngine::add({omega_a, omega_b, omega_c}))
        }
    );
    auto val_ref = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
        SymEngine::matrix_add({
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Bc, val_Dab}),
                    SymEngine::matrix_mul({minus_identity, val_Dab, val_Bc})
                }),
                val_freq_denom
            }),
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Bb, val_Dac}),
                    SymEngine::matrix_mul({minus_identity, val_Dac, val_Bb})
                }),
                val_freq_denom
            }),
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Ba, val_Dbc}),
                    SymEngine::matrix_mul({minus_identity, val_Dbc, val_Ba})
                }),
                val_freq_denom
            })
        })
    )->get_values();
    auto val_result = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
        oper_eval.apply(D->diff(a)->diff(b)->diff(c))
    )->get_values();
    REQUIRE(val_ref.size()==val_result.size());
    //FIXME: `SymEngine::expand` is necessary for the equality comparison, is it a bug for SymEngine?
    for (std::size_t i=0; i<val_ref.size(); ++i) REQUIRE(SymEngine::eq(
        *SymEngine::expand(val_ref[i]), *SymEngine::expand(val_result[i])
    ));
    // D^{aab} = 4*[Va^{a}, D^{ab}]\circ[(2*omega_a+omega_b-(E_m-E_n))^{-1}]
    //         + 2*[Vb^{b}, D^{aa}]\circ[(2*omega_a+omega_b-(E_m-E_n))^{-1}]
    auto val_Daa = oper_eval.apply(D->diff(a)->diff(a));
    auto val_freq_aba = SymEngine::immutable_dense_matrix(
        2, 2,
        {
            SymEngine::div(
                SymEngine::integer(4), SymEngine::add({omega_a, omega_a, omega_b})
            ),
            SymEngine::div(
                SymEngine::integer(4),
                SymEngine::sub(SymEngine::add({omega_a, omega_a, omega_b}), E_01)
            ),
            SymEngine::div(
                SymEngine::integer(4),
                SymEngine::sub(SymEngine::add({omega_a, omega_a, omega_b}), E_10)
            ),
            SymEngine::div(
                SymEngine::integer(4), SymEngine::add({omega_a, omega_a, omega_b})
            )
        }
    );
    auto val_freq_aab = SymEngine::immutable_dense_matrix(
        2, 2,
        {
            SymEngine::div(SymEngine::two, SymEngine::add({omega_a, omega_a, omega_b})),
            SymEngine::div(
                SymEngine::two,
                SymEngine::sub(SymEngine::add({omega_a, omega_a, omega_b}), E_01)
            ),
            SymEngine::div(
                SymEngine::two,
                SymEngine::sub(SymEngine::add({omega_a, omega_a, omega_b}), E_10)
            ),
            SymEngine::div(SymEngine::two, SymEngine::add({omega_a, omega_a, omega_b}))
        }
    );
    val_ref = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
        SymEngine::matrix_add({
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Ba, val_Dab}),
                    SymEngine::matrix_mul({minus_identity, val_Dab, val_Ba})
                }),
                val_freq_aba
            }),
            SymEngine::hadamard_product({
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({val_Bb, val_Daa}),
                    SymEngine::matrix_mul({minus_identity, val_Daa, val_Bb})
                }),
                val_freq_aab
            })
        })
    )->get_values();
    val_result = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
        oper_eval.apply(D->diff(a)->diff(a)->diff(b))
    )->get_values();
    REQUIRE(val_ref.size()==val_result.size());
    //FIXME: `SymEngine::expand` is necessary for the equality comparison, is it a bug for SymEngine?
    for (std::size_t i=0; i<val_ref.size(); ++i) REQUIRE(SymEngine::eq(
        *SymEngine::expand(val_ref[i]), *SymEngine::expand(val_result[i])
    ));

    // Test evaluator for functions
    auto fun_eval = TwoLevelFunction(
        std::make_pair(H0, val_H0),
        std::map<SymEngine::RCP<const OneElecOperator>,
                 SymEngine::RCP<const SymEngine::MatrixExpr>,
                 SymEngine::RCPBasicKeyLess>({{Va, val_Ba}, {Vb, val_Bb}, {Vc, val_Bc}}),
        std::make_pair(D, val_D)
    );

    REQUIRE(SymEngine::eq(
        *fun_eval.apply(E),
        *SymEngine::add({
            E0,
            SymEngine::mul(a, Va_00),
            SymEngine::mul(b, Vb_00),
            SymEngine::mul(c, Vc_00)
        })
    ));
    // E^{a} = tr(D*Va^{a}) + tr(D^{a}*Va) + tr(D^{a}*Vb) + tr(D^{a}*Vc)
    REQUIRE(SymEngine::eq(
        *fun_eval.apply(E->diff(a)),
        *SymEngine::add({
            Va_00,
            SymEngine::mul(
                a,
                SymEngine::add(
                    SymEngine::div(
                        SymEngine::mul({SymEngine::minus_one, Va_01, Va_10}),
                        SymEngine::sub(omega_a, E_01)
                    ),
                    SymEngine::div(
                        SymEngine::mul(Va_10, Va_01),
                        SymEngine::sub(omega_a, E_10)
                    )
                )
            ),
            SymEngine::mul(
                b,
                SymEngine::add(
                    SymEngine::div(
                        SymEngine::mul({SymEngine::minus_one, Va_01, Vb_10}),
                        SymEngine::sub(omega_a, E_01)
                    ),
                    SymEngine::div(
                        SymEngine::mul(Va_10, Vb_01), SymEngine::sub(omega_a, E_10)
                    )
                )
            ),
            SymEngine::mul(
                c,
                SymEngine::add(
                    SymEngine::div(
                        SymEngine::mul({SymEngine::minus_one, Va_01, Vc_10}),
                        SymEngine::sub(omega_a, E_01)
                    ),
                    SymEngine::div(
                        SymEngine::mul(Va_10, Vc_01), SymEngine::sub(omega_a, E_10)
                    )
                )
            )
        })
    ));
    // E^{ab} = tr(D^{b}*Va^{a}) + tr(D^{a}*Vb^{b}) + tr(D^{ab}*H0)
    //        + tr(D^{ab}*Va) + tr(D^{ab}*Vb) + tr(D^{ab}*Vc)
    REQUIRE(SymEngine::eq(
        *fun_eval.apply(E->diff(a)->diff(b)),
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({val_Db, val_Ba})),
            SymEngine::trace(SymEngine::matrix_mul({val_Da, val_Bb})),
            SymEngine::trace(SymEngine::matrix_mul({val_Dab, val_H0})),
            SymEngine::trace(SymEngine::matrix_mul({a, val_Dab, val_Ba})),
            SymEngine::trace(SymEngine::matrix_mul({b, val_Dab, val_Bb})),
            SymEngine::trace(SymEngine::matrix_mul({c, val_Dab, val_Bc}))
        })
    ));
    // E^{abc} = tr(D^{bc}*Va^{a}) + tr(D^{ac}*Vb^{b}) + tr(D^{ab}*Vc^{c})
    //         + tr(D^{abc}*H0) + tr(D^{abc}*Va) + tr(D^{abc}*Vb) + tr(D^{abc}*Vc)
    auto val_Dabc = oper_eval.apply(D->diff(a)->diff(b)->diff(c));
    REQUIRE(SymEngine::eq(
        *fun_eval.apply(E->diff(a)->diff(b)->diff(c)),
        *SymEngine::add({
            SymEngine::trace(SymEngine::matrix_mul({val_Dbc, val_Ba})),
            SymEngine::trace(SymEngine::matrix_mul({val_Dac, val_Bb})),
            SymEngine::trace(SymEngine::matrix_mul({val_Dab, val_Bc})),
            SymEngine::trace(SymEngine::matrix_mul({val_Dabc, val_H0})),
            SymEngine::trace(SymEngine::matrix_mul({a, val_Dabc, val_Ba})),
            SymEngine::trace(SymEngine::matrix_mul({b, val_Dabc, val_Bb})),
            SymEngine::trace(SymEngine::matrix_mul({c, val_Dabc, val_Bc}))
        })
    ));

    //How to test mixed_states_density()? Will V's be different?
}
