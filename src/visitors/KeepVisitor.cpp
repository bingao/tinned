#include <utility>

#include <symengine/number.h>
#include <symengine/pow.h>
#include <symengine/symengine_exception.h>

#include "Tinned/operators/PerturbedParameter.hpp"
#include "Tinned/operators/ConjugateTranspose.hpp"
#include "Tinned/operators/OneElecDensity.hpp"
#include "Tinned/operators/OneElecOperator.hpp"
#include "Tinned/operators/TwoElecEnergy.hpp"
#include "Tinned/operators/TwoElecOperator.hpp"
#include "Tinned/operators/CompositeFunction.hpp"
#include "Tinned/operators/ExchCorrEnergy.hpp"
#include "Tinned/operators/ExchCorrPotential.hpp"
#include "Tinned/operators/NonElecFunction.hpp"
#include "Tinned/operators/TemporumOperator.hpp"
#include "Tinned/operators/TemporumOverlap.hpp"
#include "Tinned/operators/AdjointMap.hpp"
#include "Tinned/operators/ClusterConjHamiltonian.hpp"

#include "Tinned/visitors/KeepVisitor.hpp"
#include "Tinned/visitors/VisitorUtilities.hpp"

namespace Tinned
{
    void KeepVisitor::bvisit(const SymEngine::Add& x)
    {
        // If `Add` will not be kept as whole, we then check its coefficient
        // and pairs
        if (condition_(x)) {
            // First we check if the coefficient will be kept
            SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
            if (condition_(*coef)) coef = SymEngine::zero;
            // Next we check each pair (`Basic` and `Number`) in the dictionary
            // of `Add`
            SymEngine::umap_basic_num d;
            for (const auto& p: x.get_dict()) {
                if (condition_(
                    *SymEngine::Add::from_dict(SymEngine::zero, {{p.first, p.second}})
                )) {
                    if (condition_(*p.second)) {
                        auto new_key = apply(p.first);
                        if (new_key.is_null()) continue;
                        SymEngine::Add::coef_dict_add_term(
                            SymEngine::outArg(coef), d, p.second, new_key
                        );
                    }
                    // `Number` will be kept so the the pair will be kept as a
                    // whole
                    else {
                        SymEngine::Add::coef_dict_add_term(
                            SymEngine::outArg(coef), d, p.second, p.first
                        );
                    }
                }
                // This pair will be kept as a whole
                else {
                    SymEngine::Add::coef_dict_add_term(
                        SymEngine::outArg(coef), d, p.second, p.first
                    );
                }
            }
            if (coef->is_zero() && d.empty()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                result_ = SymEngine::Add::from_dict(coef, std::move(d));
            }
        }
        // `Add` will be kept as a whole
        else {
            result_ = x.rcp_from_this();
        }
    }

    // This function is similar to that for MatrixMul, so see the detailed
    // explanation for that function
    void KeepVisitor::bvisit(const SymEngine::Mul& x)
    {
        if (condition_(x)) {
            SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
            // `Mul` will be kept as a whole if the coefficient is kept
            if (!condition_(*coef)) {
                result_ = x.rcp_from_this();
                return;
            }
            // Indicates if there is factor(s) kept
            bool factors_kept = false;
            // We check each pair (`Basic` and `Basic`) in the dictionary
            // of `Mul`
            SymEngine::map_basic_basic d;
            for (const auto& p : x.get_dict()) {
                // First check the whole factor, i.e. key^value
                auto factor = SymEngine::make_rcp<SymEngine::Pow>(p.first, p.second);
                if (condition_(*factor)) {
                    auto new_key = apply(p.first);
                    if (new_key.is_null()) {
                        // The key does not match any given symbols, but we
                        // save this pair in case that there will be factor(s)
                        // kept
                        SymEngine::Mul::dict_add_term_new(
                            SymEngine::outArg(coef), d, p.second, p.first
                        );
                    }
                    else {
                        // `Mul` will be kept as a whole if this key matches
                        // any given symbols and is kept
                        if (SymEngine::eq(*p.first, *new_key)) {
                            result_ = x.rcp_from_this();
                            return;
                        }
                        // (Ak+Ar)^a where Ak matches a given symbol and kept,
                        // not Ar. We simply use Newton's generalized binomial
                        // theorem
                        // (https://en.wikipedia.org/wiki/Binomial_theorem#Newton's_generalized_binomial_theorem),
                        // and save Ar^a.
                        else {
                            SymEngine::Mul::dict_add_term_new(
                                SymEngine::outArg(coef),
                                d,
                                p.second,
                                SymEngine::sub(p.first, new_key)
                            );
                            factors_kept = true;
                        }
                    }
                }
                // `Mul` will be kept as a whole if this factor matches any
                // given symbols and is kept
                else {
                    result_ = x.rcp_from_this();
                    return;
                }
            }
            // When there are factors partially kept, the result can be
            // computed as
            // c*(A^a)*(B^b)*...*(R^r)*(S^s)*... - c*(Ar^a)*(Br^b)*...*(R^r)*(S^s)*...,
            // where Ar, Br, Cr, ... are parts that are removed, R, S, T, ...
            // are those without kept parts.
            if (factors_kept) {
                result_ = SymEngine::sub(
                    x.rcp_from_this(),
                    SymEngine::Mul::from_dict(coef, std::move(d))
                );
            }
            // `Mul` will be removed since all its factors are null after
            // removal
            else {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
        }
        // `Mul` will be kept as a whole
        else {
            result_ = x.rcp_from_this();
        }
    }

    void KeepVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        // We don't allow for the removal of derivative symbols, but check only
        // if the `NonElecFunction` (or its derivatives) will be removed as a
        // whole
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            remove_if_symbol_like(SymEngine::down_cast<const NonElecFunction&>(x));
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            keep_if_a_function(
                op,
                std::bind(&construct_2el_energy, std::placeholders::_1),
                op.get_2el_operator(),
                op.get_outer_state()
            );
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            keep_if_a_function(
                op,
                std::bind(
                    &construct_composite_function,
                    std::placeholders::_1,
                    op.get_name(),
                    op.get_order()
                ),
                op.get_inner()
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            keep_if_a_function(
                op,
                std::bind(
                    &construct_xc_energy,
                    std::placeholders::_1,
                    op.get_name(),
                    op.get_state(),
                    op.get_overlap_distribution(),
                    op.get_weight()
                ),
                op.get_energy()
            );
        }
        else {
            throw SymEngine::NotImplementedError(
                "KeepVisitor::bvisit() not implemented for FunctionSymbol " + x.__str__()
            );
        }
    }

    void KeepVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const PerturbedParameter>(x)) {
            remove_if_symbol_like(SymEngine::down_cast<const PerturbedParameter&>(x));
        }
        else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
            keep_if_a_function(
                op,
                std::bind(&construct_conjugate_transpose, std::placeholders::_1),
                op.get_arg()
            );
        }
        else if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            remove_if_symbol_like(SymEngine::down_cast<const OneElecDensity&>(x));
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            remove_if_symbol_like(SymEngine::down_cast<const OneElecOperator&>(x));
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            keep_if_a_function(
                op,
                std::bind(
                    &construct_2el_operator,
                    std::placeholders::_1,
                    op.get_name(),
                    op.get_dependencies(),
                    op.get_derivatives()
                ),
                op.get_state()
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            keep_if_a_function(
                op,
                std::bind(
                    &construct_xc_potential,
                    std::placeholders::_1,
                    op.get_name(),
                    op.get_state(),
                    op.get_overlap_distribution(),
                    op.get_weight()
                ),
                op.get_potential()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            keep_if_a_function(
                op,
                std::bind(&construct_dt_operator, std::placeholders::_1, op.get_type()),
                op.get_target()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            remove_if_symbol_like(SymEngine::down_cast<const TemporumOverlap&>(x));
        }
        else if (SymEngine::is_a_sub<const AdjointMap>(x)) {
            auto& op = SymEngine::down_cast<const AdjointMap&>(x);
            keep_if_a_function(
                op,
                std::bind(&construct_adjoint_map, std::placeholders::_1),
                op.get_x(),
                op.get_y()
            );
        }
        else if (SymEngine::is_a_sub<const ClusterConjHamiltonian>(x)) {
            auto& op = SymEngine::down_cast<const ClusterConjHamiltonian&>(x);
            keep_if_a_function(
                op,
                std::bind(&construct_cc_hamiltonian, std::placeholders::_1),
                op.get_cluster_operator(),
                op.get_hamiltonian()
            );
        }
        else {
            throw SymEngine::NotImplementedError(
                "KeepVisitor::bvisit() not implemented for MatrixSymbol " + x.__str__()
            );
        }
    }

    void KeepVisitor::bvisit(const SymEngine::Trace& x)
    {
        keep_if_a_function(
            x,
            std::bind(&construct_trace, std::placeholders::_1),
            x.get_args()[0]
        );
    }

    void KeepVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        keep_if_a_function(
            x,
            std::bind(&construct_conjugate_matrix, std::placeholders::_1),
            x.get_arg()
        );
    }

    void KeepVisitor::bvisit(const SymEngine::Transpose& x)
    {
        keep_if_a_function(
            x,
            std::bind(&construct_transpose, std::placeholders::_1),
            x.get_arg()
        );
    }

    void KeepVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // If `MatrixAdd` will not be kept as whole, we then check if its
        // arguments will be kept
        if (condition_(x)) {
            SymEngine::vec_basic terms;
            for (auto& arg: x.get_args()) {
                auto new_arg = apply(arg);
                if (!new_arg.is_null()) terms.push_back(new_arg);
            }
            if (terms.empty()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                result_ = SymEngine::matrix_add(terms);
            }
        }
        // `MatrixAdd` will be kept as a whole
        else {
            result_ = x.rcp_from_this();
        }
    }

    void KeepVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // If `MatrixMul` will not be kept as whole, we then check if its
        // factors will be kept
        if (condition_(x)) {
            // `factors` will be used to construct the product that will be
            // removed. The first factor is -1 for substraction, see
            // information below.
            auto factors = SymEngine::vec_basic({SymEngine::minus_one});
            // Indicates if there is factor(s) kept
            bool factors_kept = false;
            for (auto& arg: x.get_args()) {
                auto new_arg = apply(arg);
                if (new_arg.is_null()) {
                    // This factor does not match any given symbols, but we
                    // save it in case that there will be other factor(s) kept
                    factors.push_back(arg);
                }
                else {
                    // `MatrixMul` will be kept as a whole if this factor
                    // matches a given symbol and is kept
                    if (SymEngine::eq(*arg, *new_arg)) {
                        result_ = x.rcp_from_this();
                        return;
                    }
                    else {
                        if (SymEngine::is_a_sub<const SymEngine::MatrixExpr>(*arg)) {
                            // Suppose `MatrixMul` is A*B*C*... = (Ak+Ar)*B*C*...,
                            // where Ak will be kept and Ar will be removed. The
                            // result after removal will be Ak*B*C*... . By saving
                            // Ar = A-Ak, the result can also be computed as
                            // A*B*C*... - Ar*B*C*... .
                            factors.push_back(SymEngine::matrix_add({
                                arg,
                                SymEngine::matrix_mul({SymEngine::minus_one, new_arg})
                            }));
                        }
                        // `arg` is a scalar
                        else {
                            factors.push_back(SymEngine::sub(arg, new_arg));
                        }
                        factors_kept = true;
                    }
                }
            }
            // As aforementioned, when there are factors partially kept, the
            // result can be computed as A*B*C*...*R*S*T*... - Ar*Br*Cr*...*R*S*T*...,
            // where Ar, Br, Cr, ... are parts that are removed, R, S, T, ...
            // are those without kept parts.
            if (factors_kept) {
                result_ = SymEngine::matrix_add({
                    x.rcp_from_this(), SymEngine::matrix_mul(factors)
                });
            }
            // `MatrixMul` will be removed since all its factors are null after
            // removal
            else {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
        }
        // `MatrixMul` will be kept as a whole
        else {
            result_ = x.rcp_from_this();
        }
    }
}
