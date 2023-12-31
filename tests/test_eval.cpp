#define CATCH_CONFIG_MAIN

#include <cstddef>
#include <map>
#include <string>

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
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"
#include "Tinned/FunctionEvaluator.hpp"

using namespace Tinned;

// Make a mock non-electron like function
inline SymEngine::RCP<const SymEngine::Basic> mock_nonel_function(
    const PertDependency& dependencies,
    const SymEngine::multiset_basic& derivatives
)
{
    SymEngine::vec_basic terms;
    for (const auto& dep: dependencies)
        if (dep.second>0)
            terms.push_back(SymEngine::pow(dep.first, SymEngine::integer(dep.second)));
    auto result = SymEngine::mul(terms);
    if (!derivatives.empty())
        for (const auto& var: derivatives) result = SymEngine::sdiff(result, var);
    return result;
}

// Make a mock 2x2 symmetric one-electron spin-orbital density matrix
inline SymEngine::RCP<const SymEngine::ImmutableDenseMatrix> mock_1el_density(
    const SymEngine::vec_basic& perturbations,
    const SymEngine::multiset_basic& derivatives
)
{
    SymEngine::vec_basic container;
    SymEngine::vec_basic terms;
    for (const auto& pert: perturbations)
        terms.push_back(SymEngine::pow(pert, SymEngine::minus_one));
    container.push_back(SymEngine::mul(terms));
    if (!derivatives.empty())
        for (const auto& var: derivatives)
            container.back() = SymEngine::sdiff(container.back(), var);
    container.push_back(SymEngine::mul(SymEngine::integer(3), container.front()));
    container.push_back(container.back());
    container.push_back(SymEngine::mul(SymEngine::integer(4), container.front()));
    return SymEngine::make_rcp<const SymEngine::ImmutableDenseMatrix>(2, 2, container);
}

// Make a mock 2x2 symmetric one-electron matrix
inline SymEngine::RCP<const SymEngine::ImmutableDenseMatrix> mock_1el_operator(
    const PertDependency& dependencies,
    const SymEngine::multiset_basic& derivatives
)
{
    SymEngine::vec_basic container;
    SymEngine::vec_basic terms;
    for (const auto& dep: dependencies)
        if (dep.second>0)
            terms.push_back(SymEngine::pow(dep.first, SymEngine::integer(dep.second)));
    container.push_back(SymEngine::mul(terms));
    if (!derivatives.empty())
        for (const auto& var: derivatives)
            container.back() = SymEngine::sdiff(container.back(), var);
    container.push_back(SymEngine::mul(SymEngine::integer(3), container.front()));
    container.push_back(container.back());
    container.push_back(SymEngine::mul(SymEngine::integer(4), container.front()));
    return SymEngine::make_rcp<const SymEngine::ImmutableDenseMatrix>(2, 2, container);
}

// Make mock two-electron integrals
inline std::map<std::string, SymEngine::RCP<const SymEngine::Basic>> mock_2el_integral(
    const PertDependency& dependencies,
    const SymEngine::multiset_basic& derivatives
)
{
    std::map<std::string, SymEngine::RCP<const SymEngine::Basic>> result;
    SymEngine::vec_basic terms;
    for (const auto& dep: dependencies)
        if (dep.second>0)
            terms.push_back(SymEngine::pow(dep.first, SymEngine::integer(dep.second)));
    auto eri = SymEngine::mul(terms);
    if (!derivatives.empty())
        for (const auto& var: derivatives) eri = SymEngine::sdiff(eri, var);
    // <11|11>
    result.emplace(std::string("1111"), eri);
    // <11|12>, <11|21>, <12|11>, <21|11>
    result.emplace(std::string("1112"), SymEngine::mul(SymEngine::two, eri));
    result.emplace(std::string("1121"), result[std::string("1112")]);
    result.emplace(std::string("1211"), result[std::string("1112")]);
    result.emplace(std::string("2111"), result[std::string("1112")]);
    // <11|22>, <22|11>, <21|12>, <12|21>
    result.emplace(std::string("1122"), SymEngine::mul(SymEngine::integer(3), eri));
    result.emplace(std::string("2211"), result[std::string("1122")]);
    result.emplace(std::string("2112"), result[std::string("1122")]);
    result.emplace(std::string("1221"), result[std::string("1122")]);
    // <12|12>, <21|21>
    result.emplace(std::string("1212"), SymEngine::mul(SymEngine::integer(4), eri));
    result.emplace(std::string("2121"), result[std::string("1212")]);
    // <12|22>, <21|22>, <22|12>, <22|21>
    result.emplace(std::string("1222"), SymEngine::mul(SymEngine::integer(5), eri));
    result.emplace(std::string("2122"), result[std::string("1222")]);
    result.emplace(std::string("2212"), result[std::string("1222")]);
    result.emplace(std::string("2221"), result[std::string("1222")]);
    // <22|22>
    result.emplace(std::string("2222"), SymEngine::mul(SymEngine::integer(6), eri));
    return result;
}

// Make a mock two-electron operator
inline SymEngine::RCP<const SymEngine::ImmutableDenseMatrix> mock_2el_operator(
    const PertDependency& dependencies,
    const SymEngine::multiset_basic& derivatives,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& state
)
{
    auto eri = mock_2el_integral(dependencies, derivatives);
    SymEngine::vec_basic container;
    SymEngine::vec_basic terms;
    terms.push_back(SymEngine::mul(state->get(0, 0), eri[std::string("1111")]));
    terms.push_back(SymEngine::mul(state->get(0, 1), eri[std::string("1121")]));
    terms.push_back(SymEngine::mul(state->get(1, 0), eri[std::string("1112")]));
    terms.push_back(SymEngine::mul(state->get(1, 1), eri[std::string("1122")]));
    container.push_back(SymEngine::add(terms));
    terms.clear();
    terms.push_back(SymEngine::mul(state->get(0, 0), eri[std::string("1211")]));
    terms.push_back(SymEngine::mul(state->get(0, 1), eri[std::string("1221")]));
    terms.push_back(SymEngine::mul(state->get(1, 0), eri[std::string("1212")]));
    terms.push_back(SymEngine::mul(state->get(1, 1), eri[std::string("1222")]));
    container.push_back(SymEngine::add(terms));
    terms.clear();
    terms.push_back(SymEngine::mul(state->get(0, 0), eri[std::string("2111")]));
    terms.push_back(SymEngine::mul(state->get(0, 1), eri[std::string("2121")]));
    terms.push_back(SymEngine::mul(state->get(1, 0), eri[std::string("2112")]));
    terms.push_back(SymEngine::mul(state->get(1, 1), eri[std::string("2122")]));
    container.push_back(SymEngine::add(terms));
    terms.clear();
    terms.push_back(SymEngine::mul(state->get(0, 0), eri[std::string("2211")]));
    terms.push_back(SymEngine::mul(state->get(0, 1), eri[std::string("2221")]));
    terms.push_back(SymEngine::mul(state->get(1, 0), eri[std::string("2212")]));
    terms.push_back(SymEngine::mul(state->get(1, 1), eri[std::string("2222")]));
    container.push_back(SymEngine::add(terms));
    return SymEngine::make_rcp<const SymEngine::ImmutableDenseMatrix>(2, 2, container);
}

// Make a mock two-electron energy
inline SymEngine::RCP<const SymEngine::Basic> mock_2el_energy(
    const PertDependency& dependencies,
    const SymEngine::multiset_basic& derivatives,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& inner,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& outer
)
{
    auto G = mock_2el_operator(dependencies, derivatives, inner);
std::cout << "G = " << stringify(G) << "\n";
    return SymEngine::trace(SymEngine::matrix_mul(SymEngine::vec_basic({G, outer})));
}

// Make a mock generalized density vector
inline SymEngine::RCP<const SymEngine::Basic> mock_exc_density(
    const SymEngine::RCP<const SymEngine::Trace>& x,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& state,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& Omega
)
{
    auto arg = x->get_args()[0];
    if (SymEngine::is_a_sub<const SymEngine::MatrixMul>(*arg)) {
        auto op = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(arg);
        auto factors = op->get_args();
        if (factors.size() != 2) {
            throw SymEngine::SymEngineException(
                "Invalid number of arguments from generalized density vectors "
                + stringify(x)
            );
        }
        SymEngine::vec_basic terms;
        for (const auto& factor: factors) {
            SymEngine::multiset_basic derivatives;
            if (SymEngine::is_a_sub<const OneElecDensity>(*factor)) {
                terms.push_back(state);
                derivatives = SymEngine::rcp_dynamic_cast<const OneElecDensity>(factor)->get_derivatives();
            }
            else if (SymEngine::is_a_sub<const OneElecOperator>(*factor)) {
                terms.push_back(Omega);
                derivatives = SymEngine::rcp_dynamic_cast<const OneElecOperator>(factor)->get_derivatives();
            }
            else {
                throw SymEngine::SymEngineException(
                    "Invalid argument "
                    + stringify(factor)
                    + " from generalized density vectors "
                    + stringify(x)
                );
            }
            if (!derivatives.empty())
                for (const auto& var: derivatives)
                    terms.back() = SymEngine::sdiff(terms.back(), var);
        }
        return SymEngine::trace(SymEngine::matrix_mul(terms));
    }
    else {
        throw SymEngine::SymEngineException(
            "Invalid generalized density vectors " + stringify(x)
        );
    }
}

// Make a mock exchange-correlation energy
inline SymEngine::RCP<const SymEngine::Basic> mock_xc_energy(
    const ExchCorrEnergy& x,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& state,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& Omega,
    const SymEngine::RCP<const SymEngine::Basic>& weight
)
{
    auto density = SymEngine::trace(SymEngine::matrix_mul({Omega, state}));
    SymEngine::vec_basic energy_terms;
    for (const auto& energy_term: x.get_energy_terms()) {
        SymEngine::vec_basic factors;
        for (const auto& arg: energy_term->get_args()) {
            if (SymEngine::is_a_sub<const SymEngine::Number>(*arg)) {
                factors.push_back(arg);
            }
            // Grid weight
            else if (SymEngine::is_a_sub<const NonElecFunction>(*arg)) {
                factors.push_back(weight);
                auto op = SymEngine::rcp_dynamic_cast<const NonElecFunction>(arg);
                auto derivatives = op->get_derivatives();
                if (!derivatives.empty())
                    for (const auto& var: derivatives)
                        factors.back() = SymEngine::sdiff(factors.back(), var);
            }
            // XC functional derivative
            else if (SymEngine::is_a_sub<const CompositeFunction>(*arg)) {
                auto op = SymEngine::rcp_dynamic_cast<const CompositeFunction>(arg);
                factors.push_back(
                    SymEngine::pow(density, SymEngine::integer(-op->get_order()-1))
                );
            }
            // Generalized density vector, sum of generalized density
            // vectors, or power of (sum of) generalized density
            // vector(s)
            else if (SymEngine::is_a_sub<const SymEngine::Trace>(*arg)) {
                factors.push_back(mock_exc_density(
                    SymEngine::rcp_dynamic_cast<const SymEngine::Trace>(arg),
                    state,
                    Omega
                ));
            }
            else if (SymEngine::is_a_sub<const SymEngine::Add>(*arg)) {
                auto op = SymEngine::rcp_dynamic_cast<const SymEngine::Add>(arg);
                SymEngine::vec_basic sum_terms;
                for (const auto& term: op->get_args()) {
                    if (SymEngine::is_a_sub<const SymEngine::Trace>(*term)) {
                        sum_terms.push_back(mock_exc_density(
                            SymEngine::rcp_dynamic_cast<const SymEngine::Trace>(term),
                            state,
                            Omega
                        ));
                    }
                    else {
                        throw SymEngine::SymEngineException(
                            "Invalid generalized density vector "
                            + stringify(term)
                            + " from the XC energy density contraction "
                            + stringify(energy_term)
                        );
                    }
                }
                factors.push_back(SymEngine::add(sum_terms));
            }
            else if (SymEngine::is_a_sub<const SymEngine::Pow>(*arg)) {
                auto op = SymEngine::rcp_dynamic_cast<const SymEngine::Pow>(arg);
                auto base = op->get_base();
                auto exponent = op->get_exp();
                if (!SymEngine::is_a_sub<const SymEngine::Trace>(*base)) {
                    throw SymEngine::SymEngineException(
                        "Invalid base of generalized density vector "
                        + stringify(arg)
                        + " from the XC energy density contraction "
                        + stringify(energy_term)
                    );
                }
                if (!SymEngine::is_a_sub<const SymEngine::Number>(*exponent)) {
                    throw SymEngine::SymEngineException(
                        "Invalid exponent of generalized density vector "
                        + stringify(arg)
                        + " from the XC energy density contraction "
                        + stringify(energy_term)
                    );
                }
                factors.push_back(SymEngine::pow(
                    mock_exc_density(
                        SymEngine::rcp_dynamic_cast<const SymEngine::Trace>(base),
                        state,
                        Omega
                    ),
                    exponent
                ));
            }
            else {
                throw SymEngine::SymEngineException(
                    "Invalid type from the XC energy density contraction "
                    + stringify(energy_term)
                );
            }
        }
        energy_terms.push_back(SymEngine::mul(factors));
    }
    return SymEngine::add(energy_terms);
}

// Make a mock exchange-correlation potential
inline SymEngine::RCP<const SymEngine::Basic> mock_xc_potential(
    const ExchCorrPotential& x,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& state,
    const SymEngine::RCP<const SymEngine::ImmutableDenseMatrix>& Omega,
    const SymEngine::RCP<const SymEngine::Basic>& weight
)
{
    auto density = SymEngine::trace(SymEngine::matrix_mul({Omega, state}));
    SymEngine::vec_basic potential_terms;
    for (const auto& potential_term: x.get_potential_terms()) {
        SymEngine::vec_basic factors;
        for (const auto& arg: potential_term->get_args()) {
            // Overlap distribution
            if (SymEngine::is_a_sub<const OneElecOperator>(*arg)) {
                factors.push_back(Omega);
                auto op = SymEngine::rcp_dynamic_cast<const OneElecOperator>(arg);
                auto derivatives = op->get_derivatives();
                if (!derivatives.empty())
                    for (const auto& var: derivatives)
                        factors.back() = SymEngine::sdiff(factors.back(), var);
            }
            else if (SymEngine::is_a_sub<const SymEngine::Number>(*arg)) {
                factors.push_back(arg);
            }
            // Grid weight
            else if (SymEngine::is_a_sub<const NonElecFunction>(*arg)) {
                factors.push_back(weight);
                auto op = SymEngine::rcp_dynamic_cast<const NonElecFunction>(arg);
                auto derivatives = op->get_derivatives();
                if (!derivatives.empty())
                    for (const auto& var: derivatives)
                        factors.back() = SymEngine::sdiff(factors.back(), var);
            }
            // XC functional derivative
            else if (SymEngine::is_a_sub<const CompositeFunction>(*arg)) {
                auto op = SymEngine::rcp_dynamic_cast<const CompositeFunction>(arg);
                factors.push_back(
                    SymEngine::pow(density, SymEngine::integer(-op->get_order()-1))
                );
            }
            // Generalized density vector, sum of generalized density
            // vectors, or power of (sum of) generalized density
            // vector(s)
            else if (SymEngine::is_a_sub<const SymEngine::Trace>(*arg)) {
                factors.push_back(mock_exc_density(
                    SymEngine::rcp_dynamic_cast<const SymEngine::Trace>(arg),
                    state,
                    Omega
                ));
            }
            else if (SymEngine::is_a_sub<const SymEngine::Add>(*arg)) {
                auto op = SymEngine::rcp_dynamic_cast<const SymEngine::Add>(arg);
                SymEngine::vec_basic sum_terms;
                for (const auto& term: op->get_args()) {
                    if (SymEngine::is_a_sub<const SymEngine::Trace>(*term)) {
                        sum_terms.push_back(mock_exc_density(
                            SymEngine::rcp_dynamic_cast<const SymEngine::Trace>(term),
                            state,
                            Omega
                        ));
                    }
                    else {
                        throw SymEngine::SymEngineException(
                            "Invalid generalized density vector "
                            + stringify(term)
                            + " from the XC potential density contraction "
                            + stringify(potential_term)
                        );
                    }
                }
                factors.push_back(SymEngine::add(sum_terms));
            }
            else if (SymEngine::is_a_sub<const SymEngine::Pow>(*arg)) {
                auto op = SymEngine::rcp_dynamic_cast<const SymEngine::Pow>(arg);
                auto base = op->get_base();
                auto exponent = op->get_exp();
                if (!SymEngine::is_a_sub<const SymEngine::Trace>(*base)) {
                    throw SymEngine::SymEngineException(
                        "Invalid base of generalized density vector "
                        + stringify(arg)
                        + " from the XC potential density contraction "
                        + stringify(potential_term)
                    );
                }
                if (!SymEngine::is_a_sub<const SymEngine::Number>(*exponent)) {
                    throw SymEngine::SymEngineException(
                        "Invalid exponent of generalized density vector "
                        + stringify(arg)
                        + " from the XC potential density contraction "
                        + stringify(potential_term)
                    );
                }
                factors.push_back(SymEngine::pow(
                    mock_exc_density(
                        SymEngine::rcp_dynamic_cast<const SymEngine::Trace>(base),
                        state,
                        Omega
                    ),
                    exponent
                ));
            }
            else {
                throw SymEngine::SymEngineException(
                    "Invalid type from the XC potential density contraction "
                    + stringify(potential_term)
                );
            }
        }
        potential_terms.push_back(SymEngine::matrix_mul(factors));
    }
    return SymEngine::matrix_add(potential_terms);
}

class FunMockEvaluator: public FunctionEvaluator<SymEngine::RCP<const SymEngine::Basic>>
{
    private:
        SymEngine::vec_basic perturbations_;

    protected:
        SymEngine::RCP<const SymEngine::Basic> eval_nonel_function(
            const NonElecFunction& x
        ) override
        {
            return mock_nonel_function(x.get_dependencies(), x.get_derivatives());
        }

        SymEngine::RCP<const SymEngine::Basic> eval_2el_energy(
            const TwoElecEnergy& x
        ) override
        {
std::cout << "x = " << stringify(x) << "\n";
            auto inner = mock_1el_density(
                perturbations_,
                x.get_inner_state()->get_derivatives()
            );
            auto outer = mock_1el_density(
                perturbations_,
                x.get_outer_state()->get_derivatives()
            );
std::cout << "inner = " << stringify(inner) << "\n";
std::cout << "outer = " << stringify(outer) << "\n";
            return mock_2el_energy(
                x.get_dependencies(), x.get_derivatives(), inner, outer
            );
        }

        SymEngine::RCP<const SymEngine::Basic> eval_xc_energy(
            const ExchCorrEnergy& x
        ) override
        {
            auto state = mock_1el_density(
                perturbations_,
                x.get_state()->get_derivatives()
            );
            auto Omega = mock_1el_operator(
                x.get_overlap_distribution()->get_dependencies(),
                x.get_overlap_distribution()->get_derivatives()
            );
            auto weight = mock_nonel_function(
                x.get_weight()->get_dependencies(),
                x.get_weight()->get_derivatives()
            );
            return mock_xc_energy(x, state, Omega, weight);
        }

        SymEngine::RCP<const SymEngine::Basic> eval_trace(
            const SymEngine::RCP<const SymEngine::Basic>& scalar,
            const SymEngine::vec_basic& factors
        ) override
        {
            if (SymEngine::is_a_Number(*scalar)) {
                SymEngine::vec_basic terms;
                for (const auto& factor: factors) {
                    if (SymEngine::is_a_sub<const OneElecDensity>(*factor)) {
                        auto op = SymEngine::rcp_dynamic_cast<const OneElecDensity>(factor);
                        terms.push_back(
                            mock_1el_density(perturbations_, op->get_derivatives())
                        );
                    }
                    else if (SymEngine::is_a_sub<const OneElecOperator>(*factor)) {
                        auto op = SymEngine::rcp_dynamic_cast<const OneElecOperator>(factor);
                        terms.push_back(mock_1el_operator(
                            op->get_dependencies(), op->get_derivatives()
                        ));
                    }
                    else if (SymEngine::is_a_sub<const TwoElecOperator>(*factor)) {
                        auto op = SymEngine::rcp_dynamic_cast<const TwoElecOperator>(factor);
                        terms.push_back(mock_2el_operator(
                            op->get_dependencies(),
                            op->get_derivatives(),
                            mock_1el_density(
                                perturbations_, op->get_state()->get_derivatives()
                            )
                        ));
                    }
                    else if (SymEngine::is_a_sub<const ExchCorrPotential>(*factor)) {
                        auto op = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(factor);
                        terms.push_back(mock_xc_potential(
                            *op,
                            mock_1el_density(
                                perturbations_,
                                op->get_state()->get_derivatives()
                            ),
                            mock_1el_operator(
                                op->get_overlap_distribution()->get_dependencies(),
                                op->get_overlap_distribution()->get_derivatives()
                            ),
                            mock_nonel_function(
                                op->get_weight()->get_dependencies(),
                                op->get_weight()->get_derivatives()
                            )
                        ));
                    }
                    else if (SymEngine::is_a_sub<const TemporumOperator>(*factor)) {
                        auto op = SymEngine::rcp_dynamic_cast<const TemporumOperator>(factor);
                        auto target = op->get_target();
                        if (SymEngine::is_a_sub<const OneElecDensity>(*target)) {
                            terms.push_back(
                                mock_1el_density(perturbations_, op->get_derivatives())
                            );
                        }
                        else if (SymEngine::is_a_sub<const OneElecOperator>(*target)) {
                            auto op_target = SymEngine::rcp_dynamic_cast<const OneElecOperator>(target);
                            terms.push_back(mock_1el_operator(
                                op_target->get_dependencies(), op_target->get_derivatives()
                            ));
                        }
                        else {
                            throw SymEngine::NotImplementedError(
                                "Not implemented for TemporumOperator target "
                                + stringify(target)
                            );
                        }
                        terms.back() = SymEngine::mul(op->get_frequency(), terms.back());
                    }
                    else if (SymEngine::is_a_sub<const TemporumOverlap>(*factor)) {
                        auto op = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(factor);
                        SymEngine::vec_basic braket;
                        for (std::size_t i = 0; i<op->size(); ++i) {
                            // Notice that this is not the correct half
                            // time-differentiated overlap integrals
                            braket.push_back(SymEngine::mul(
                                op->get_frequency(i),
                                mock_1el_operator(
                                    op->get_dependencies(), op->get_derivatives(i).first
                                )
                            ));
                        }
                        terms.push_back(SymEngine::add(braket));
                    }
                    else {
                        throw SymEngine::NotImplementedError(
                            "Not implemented for MatrixSymbol " + stringify(factor)
                        );
                    }
                }
                return SymEngine::trace(SymEngine::matrix_mul(terms));
            }
            else {
                throw SymEngine::NotImplementedError(
                    "Non-number " + stringify(scalar) + " is not supported!"
                );
            }
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
        explicit FunMockEvaluator(const SymEngine::vec_basic& perturbations):
            perturbations_(perturbations) {}
        ~FunMockEvaluator() = default;
};

#include <iostream>

TEST_CASE("Test FunctionEvaluator", "[FunctionEvaluator]")
{
    auto a = make_perturbation(std::string("a"));
    auto b = make_perturbation(std::string("b"));
    auto dependencies = PertDependency({std::make_pair(a, 9), std::make_pair(b, 9)});
    auto D = make_1el_density(std::string("D"));
    auto h = make_1el_operator(std::string("h"), dependencies);
    auto weight = make_nonel_function(std::string("weight"));
    auto Omega = make_1el_operator(std::string("Omega"), dependencies);
    auto Exc = make_xc_energy(std::string("Exc"), D, Omega, weight);
    auto hnuc = make_nonel_function(std::string("hnuc"), dependencies);
    auto E = SymEngine::add(SymEngine::vec_basic({
        SymEngine::trace(SymEngine::matrix_mul(SymEngine::vec_basic({h, D}))),
        make_2el_energy(std::string("G"), D, D, dependencies),
        //Exc,
        hnuc
    }));

    auto evaluator = FunMockEvaluator(SymEngine::vec_basic({a, b}));
    std::cout << stringify(evaluator.apply(E)) << "\n";
}
