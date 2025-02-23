#include <utility>

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

#include "Tinned/visitors/RemoveVisitor.hpp"

namespace Tinned
{
    void RemoveVisitor::bvisit(const SymEngine::Basic& x)
    {
        throw SymEngine::NotImplementedError(
            "RemoveVisitor::bvisit() not implemented for "+x.__str__()
        );
    }

    void RemoveVisitor::bvisit(const SymEngine::Symbol& x)
    {
        remove_if_symbol_like(x);
    }

    void RemoveVisitor::bvisit(const SymEngine::Number& x)
    {
        remove_if_symbol_like(x);
    }

    void RemoveVisitor::bvisit(const SymEngine::Add& x)
    {
        // We first check if `Add` will be removed as a whole
        if (condition_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            // Indicate if `Add` will be kept as a whole
            bool kept = true;
            // First we check if the coefficient will be removed
            SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
            if (condition_(*coef)) {
                kept = false;
                coef = SymEngine::zero;
            }
            // Next we check each pair (`Basic` and `Number`) in the dictionary
            // of `Add`
            SymEngine::umap_basic_num d;
            for (const auto& p: x.get_dict()) {
                // Skip if this pair will be removed as a whole
                if (condition_(
                    *SymEngine::Add::from_dict(SymEngine::zero, {{p.first, p.second}})
                )) {
                    kept = false;
                    continue;
                }
                // Skip if `Basic` was removed
                auto new_key = apply(p.first);
                if (new_key.is_null()) {
                    kept = false;
                    continue;
                }
                // Skip if `Number` will be removed
                if (condition_(*p.second)) {
                    kept = false;
                    continue;
                }
                if (SymEngine::neq(*p.first, *new_key)) kept = false;
                SymEngine::Add::coef_dict_add_term(
                    SymEngine::outArg(coef), d, p.second, new_key
                );
            }
            if (kept) {
                result_ = x.rcp_from_this();
            }
            else {
                if (coef->is_zero() && d.empty()) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    result_ = SymEngine::Add::from_dict(coef, std::move(d));
                }
            }
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Mul& x)
    {
        // We first check if `Mul` will be removed as a whole
        if (condition_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            // We remove the whole `Mul` if its coefficient will be removed
            SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
            if (condition_(*coef)) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                // Indicate if `Mul` will be kept as a whole
                bool kept = true;
                // We check each pair (`Basic` and `Basic`) in the dictionary
                // of `Mul`
                SymEngine::map_basic_basic d;
                for (const auto& p : x.get_dict()) {
                    // We remove the whole `Mul` if the pair will be removed
                    auto factor = SymEngine::make_rcp<SymEngine::Pow>(p.first, p.second);
                    if (condition_(*factor)) {
                        result_ = SymEngine::RCP<const SymEngine::Basic>();
                        return;
                    }
                    // Remove the whole `Mul` if the key will be removed
                    auto new_key = apply(p.first);
                    if (new_key.is_null()) {
                        result_ = SymEngine::RCP<const SymEngine::Basic>();
                        return;
                    }
                    // The value in the pair (the exponent) is not allowed to
                    // be removed completely
                    auto new_value = apply(p.second);
                    if (new_value.is_null()) throw SymEngine::SymEngineException(
                        "RemoveVisitor::bvisit() does not allow to remove the exponent in a key-value pair of Mul."
                    );
                    if (SymEngine::neq(*p.first, *new_key) ||
                        SymEngine::neq(*p.second, *new_value)) kept = false;
                    SymEngine::Mul::dict_add_term_new(
                        SymEngine::outArg(coef), d, new_value, new_key
                    );
                }
                result_ = kept
                        ? x.rcp_from_this()
                        : SymEngine::Mul::from_dict(coef, std::move(d));
            }
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Constant& x)
    {
        remove_if_symbol_like(x);
    }

    void RemoveVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        // We don't allow for the removal of derivative symbols, but check only
        // if the `NonElecFunction` (or its derivatives) will be removed as a
        // whole
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            remove_if_symbol_like(SymEngine::down_cast<const NonElecFunction&>(x));
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            remove_if_a_function(
                op,
                std::bind(&construct_2el_energy, std::placeholders::_1),
                op.get_2el_operator(),
                op.get_outer_state()
            );
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            remove_if_a_function(
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
            remove_if_a_function(
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
                "RemoveVisitor::bvisit() not implemented for FunctionSymbol "+x.__str__()
            );
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        remove_if_symbol_like(x);
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const PerturbedParameter>(x)) {
            remove_if_symbol_like(SymEngine::down_cast<const PerturbedParameter&>(x));
        }
        else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
            remove_if_a_function(
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
            remove_if_a_function(
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
            remove_if_a_function(
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
            remove_if_a_function(
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
            remove_if_a_function(
                op,
                std::bind(&construct_adjoint_map, std::placeholders::_1),
                op.get_x(),
                op.get_y()
            );
        }
        else if (SymEngine::is_a_sub<const ClusterConjHamiltonian>(x)) {
            auto& op = SymEngine::down_cast<const ClusterConjHamiltonian&>(x);
            remove_if_a_function(
                op,
                std::bind(&construct_cc_hamiltonian, std::placeholders::_1),
                op.get_cluster_operator(),
                op.get_hamiltonian()
            );
        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveVisitor::bvisit() not implemented for MatrixSymbol "+x.__str__()
            );
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Trace& x)
    {
        remove_if_a_function(
            x,
            std::bind(&construct_trace, std::placeholders::_1),
            x.get_args()[0]
        );
    }

    void RemoveVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        remove_if_a_function(
            x,
            std::bind(&construct_conjugate_matrix, std::placeholders::_1),
            x.get_arg()
        );
    }

    void RemoveVisitor::bvisit(const SymEngine::Transpose& x)
    {
        remove_if_a_function(
            x,
            std::bind(&construct_transpose, std::placeholders::_1),
            x.get_arg()
        );
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We first check if `MatrixAdd` will be removed as a whole
        if (condition_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            // Indicate if `MatrixAdd` will be kept as a whole
            bool kept = true;
            // Next we check each argument of `MatrixAdd`
            SymEngine::vec_basic terms;
            for (auto& arg: x.get_args()) {
                auto new_arg = apply(arg);
                if (new_arg.is_null()) {
                    kept = false;
                }
                else {
                    if (SymEngine::neq(*arg, *new_arg)) kept = false;
                    terms.push_back(new_arg);
                }
            }
            if (terms.empty()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                if (kept) {
                    result_ = x.rcp_from_this();
                }
                else {
                    result_ = SymEngine::matrix_add(terms);
                }
            }
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // We first check if `MatrixMul` will be removed as a whole
        if (condition_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            // Indicate if `MatrixMul` will be kept as a whole
            bool kept = true;
            // Next we check each argument of `MatrixMul`
            SymEngine::vec_basic factors;
            for (auto& arg: x.get_args()) {
                auto new_arg = apply(arg);
                if (new_arg.is_null()) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                    return;
                }
                else {
                    if (SymEngine::neq(*arg, *new_arg)) kept = false;
                    factors.push_back(new_arg);
                }
            }
            // Probably we do not need to check if factors is empty because the
            // above loop over arguments should be at least executed once
            if (factors.empty()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                if (kept) {
                    result_ = x.rcp_from_this();
                }
                else {
                    result_ = SymEngine::matrix_mul(factors);
                }
            }
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        // `MatrixDerivative` represents derivatives of a `MatrixSymbol`
        // object, so according to rule "(3) Symbols and their derivatives are
        // different for the removal procedure", we need only check if
        // `MatrixDerivative` will be removed as a whole
        remove_if_symbol_like(x);
    }
}
