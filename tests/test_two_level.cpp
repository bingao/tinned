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
#include "Tinned/OperatorEvaluator.hpp"
#include "Tinned/FunctionEvaluator.hpp"

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

// Evaluator for different operators
class TwoLevelOperator: public OperatorEvaluator<SymEngine::RCP<const SymEngine::MatrixExpr>>
{
    protected:
        // Maximum allowed order
        unsigned int max_order_;
        // Unperturbed Hamiltonian and its value
        std::pair<SymEngine::RCP<const SymEngine::Basic>,
                  SymEngine::RCP<const SymEngine::MatrixExpr>> H0_;
        // Perturbation, and corresponding external field's operator and its value
        std::map<SymEngine::RCP<const Perturbation>,
                 std::pair<SymEngine::RCP<const SymEngine::Basic>,
                           SymEngine::RCP<const SymEngine::MatrixExpr>>,
                 SymEngine::RCPBasicKeyLess> V_;
        // Density matrix and the value of unperturbed one
        std::pair<SymEngine::RCP<const SymEngine::Basic>,
                  SymEngine::RCP<const SymEngine::MatrixExpr>> rho_;
        // All perturbations
        SymEngine::set_basic perturbations_;
        // Cached derivatives of density matrix <order, <perturbations, derivatives>>
        //FIXME: not optimal for identical perturbations
        std::map<unsigned int,
                 std::vector<std::pair<SymEngine::vec_basic,
                                       SymEngine::RCP<const SymEngine::MatrixExpr>>>>
            rho_derivatives_;

        SymEngine::vec_basic eval_pert_parameter(const PerturbedParameter& x) override
        {
            throw SymEngine::SymEngineException(
                "PerturbedParameter is not allowed: "+stringify(x)
            );
        }

        SymEngine::RCP<const SymEngine::MatrixExpr>
        eval_hermitian_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override
        {
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
            SymEngine::vec_basic values;
            for (std::size_t j=0; j<op->ncols(); ++j)
                for (std::size_t i=0; i<op->nrows(); ++i)
                    values.push_back(op->get(i, j)->conjugate());
            return SymEngine::immutable_dense_matrix(op->ncols(), op->nrows(), values);
        }

        SymEngine::RCP<const SymEngine::MatrixExpr>
        eval_1el_density(const OneElecDensity& x) override
        {
            auto derivatives = x.get_derivatives();
            if (derivatives.size()>max_order_) SymEngine::SymEngineException(
                "Invalid order: "+std::to_string(derivatives.size())+"("+std::to_string(max_order_)+")"
            );
            switch (derivatives.size()) {
                case 0:
                    return rho_.second;
                default:
                    if (derivatives.size()>rho_derivatives_.size()) {
                        for (unsigned int order=rho_derivatives_.size(); order<derivatives.size(); ++order) {
                            std::vector<std::pair<SymEngine::vec_basic,
                                                  SymEngine::RCP<const SymEngine::MatrixExpr>>>
                                rho_curr_derivatives;
                            auto pert_permutation = PertPermutation(order+1, perturbations_);
                            bool remaining = true;
                            do {
                                auto permut_derivative = pert_permutation.get_derivatives(remaining);
                                auto freq_sum = get_frequency_sum(permut_derivative);
                                // The last perturbation is the one for the external field
                                auto pert_field = permut_derivative.back();
                                // Find the lower order derivative of density matrix
                                for (const auto& derivative: rho_derivatives_[order]) {
                                    bool lower_found = true;
                                    for (std::size_t i=0; i<derivative.first.size(); ++i)
                                        if (SymEngine::neq(permut_derivative[i], derivative.first[i])) {
                                            lower_found = false;
                                            break;
                                        }
                                    if (lower_found) {
                                        auto val_rho = SymEngine::matrix_add({
                                            derivative.second
                                        });
                                        rho_curr_derivatives.push_back(permut_derivative, val_rho);
                                        break;
                                    }
                                }
                                if () throw
                            } while (remaining);

                            rho_derivatives_[order+1] = rho_curr_derivatives;
                        }
                    }
                    for (const auto& derivative: rho_derivatives_[derivatives.size()]) {
                        if (SymEngine::unified_eq(derivatives, derivative.first))
                            return derivative.second;
                    }
                    SymEngine::SymEngineException(
                        "Density matrix with invalid derivatives : "+stringify(x)
                    );
            }
        }

        SymEngine::RCP<const SymEngine::MatrixExpr>
        eval_1el_operator(const OneElecOperator& x) override
        {
            if (x.get_name())
        }

        SymEngine::RCP<const SymEngine::MatrixExpr>
        eval_2el_operator(const TwoElecOperator& x) override
        {
            throw SymEngine::SymEngineException(
                "TwoElecOperator is not allowed: "+stringify(x)
            );
        }

        SymEngine::RCP<const SymEngine::MatrixExpr>
        eval_temporum_operator(const TemporumOperator& x) override
        {
            result_ = apply(x.get_target());
            eval_oper_scale(x.get_frequency(), result_);
        }

        SymEngine::RCP<const SymEngine::MatrixExpr>
        eval_temporum_overlap(const TemporumOverlap& x) override
        {
            throw SymEngine::SymEngineException(
                "TemporumOverlap is not allowed: "+stringify(x)
            );
        }

        SymEngine::RCP<const SymEngine::MatrixExpr>
        eval_conjugate_matrix(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override
        {
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
            SymEngine::vec_basic values;
            for (std::size_t i=0; i<op->nrows(); ++i)
                for (std::size_t j=0; j<op->ncols(); ++j)
                    values.push_back(op->get(i, j)->conjugate());
            return SymEngine::immutable_dense_matrix(op->nrows(), op->ncols(), values);
        }

        SymEngine::RCP<const SymEngine::MatrixExpr>
        eval_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A) override
        {
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
            SymEngine::vec_basic values;
            for (std::size_t j=0; j<op->ncols(); ++j)
                for (std::size_t i=0; i<op->nrows(); ++i)
                    values.push_back(op->get(i, j));
            return SymEngine::immutable_dense_matrix(op->ncols(), op->nrows(), values);
        }

        void eval_oper_addition(
            SymEngine::RCP<const SymEngine::MatrixExpr>& A,
            const SymEngine::RCP<const SymEngine::MatrixExpr>& B
        ) override
        {
            auto op_A = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
            auto op_B = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(B);
            SYMENGINE_ASSERT(
                op_A->nrows()==op_B->nrows() && op_A->ncols()==op_B->ncols()
            )
            SymEngine::vec_basic values;
            for (std::size_t i=0; i<op_A->nrows(); ++i)
                for (std::size_t j=0; j<op_A->ncols(); ++j)
                    values.push_back(
                        SymEngine::add(op_A->get(i, j), op_B->get(i, j))
                    );
            A = SymEngine::immutable_dense_matrix(op_A->nrows(), op_A->ncols(), values);
        }

        SymEngine::RCP<const SymEngine::MatrixExpr> eval_oper_multiplication(
            const SymEngine::RCP<const SymEngine::MatrixExpr>& A,
            const SymEngine::RCP<const SymEngine::MatrixExpr>& B
        ) override
        {
            return SymEngine::matrix_mul({A, B});
        }

        void eval_oper_scale(
            const SymEngine::RCP<const SymEngine::Number>& scalar,
            SymEngine::RCP<const SymEngine::MatrixExpr>& A
        ) override
        {
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
            SymEngine::vec_basic values;
            for (std::size_t i=0; i<op->nrows(); ++i)
                for (std::size_t j=0; j<op->ncols(); ++j)
                    values.push_back(SymEngine::mul(scalar, op->get(i, j)));
            A = SymEngine::immutable_dense_matrix(op->nrows(), op->ncols(), values);
        }

    public:
        explicit TwoLevelOperator(
            const SymEngine::RCP<const SymEngine::MatrixExpr>& H0,
            const std::map<SymEngine::RCP<const Perturbation>,
                           SymEngine::RCP<const SymEngine::MatrixExpr>,
                           SymEngine::RCPBasicKeyLess()>& V_mapping,
            const SymEngine::RCP<const SymEngine::MatrixExpr>& rho0
        ) : max_order_(7), H0_(H0), V_mapping_(V_mapping)
        {
            // Density matrix should be idempotent and has purity one
            auto Z = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
                SymEngine::matrix_add({
                    SymEngine::matrix_mul({rho0, rho0}),
                    SymEngine::matrix_mul({SymEngine::minus_one, rho0})
                })
            );
            for (const auto& val: Z->get_values())
                if (!is_zero_quantity(val)) throw SymEngine::SymEngineException(
                    "Density matrix isn't idempotent: "+stringify(rho0)
                );
            if (!is_zero_quantity(SymEngine::add(SymEngine::minus_one, SymEngine::trace(rho0))))
                throw SymEngine::SymEngineException(
                    "Density matrix doesn't have purity one: "+stringify(rho0)
                );
            rho0_ = rho0;
        }

        ~TwoLevelOperator() = default;
}

class TwoLevelFunction: public FunctionEvaluator<SymEngine::RCP<const SymEngine::Basic>>
{
    private:
        SymEngine::vec_basic perturbations_;

    protected:
        SymEngine::RCP<const SymEngine::Basic> eval_nonel_function(
            const NonElecFunction& x
        ) override
        {
        }

        SymEngine::RCP<const SymEngine::Basic> eval_2el_energy(
            const TwoElecEnergy& x
        ) override
        {
        }

        SymEngine::RCP<const SymEngine::Basic> eval_xc_energy(
            const ExchCorrEnergy& x
        ) override
        {
        }

        SymEngine::RCP<const SymEngine::Basic> eval_trace(
            const SymEngine::RCP<const SymEngine::Basic>& scalar,
            const SymEngine::vec_basic& factors
        ) override
        {
        }

        void eval_fun_addition(
            SymEngine::RCP<const SymEngine::Basic>& f,
            const SymEngine::RCP<const SymEngine::Basic>& g
        ) override
        {
            f = SymEngine::add(f, g);
        }

        void eval_fun_scale(
            const SymEngine::RCP<const SymEngine::Number>& scalar,
            SymEngine::RCP<const SymEngine::Basic>& fun
        ) override
        {
            fun = SymEngine::mul(scalar, fun);
        }

    public:
        explicit TwoLevelFunction(const SymEngine::vec_basic& perturbations):
            perturbations_(perturbations) {}
        ~TwoLevelFunction() = default;
};

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
