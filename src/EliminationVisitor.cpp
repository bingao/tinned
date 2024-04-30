#include <utility>

#include <symengine/symengine_exception.h>

#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecEnergy.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

#include "Tinned/LagMultiplier.hpp"
#include "Tinned/StateVector.hpp"
#include "Tinned/StateOperator.hpp"
#include "Tinned/AdjointMap.hpp"
#include "Tinned/ExpAdjointHamiltonian.hpp"

#include "Tinned/EliminationVisitor.hpp"

namespace Tinned
{
    void EliminationVisitor::bvisit(const SymEngine::Basic& x)
    {
        throw SymEngine::NotImplementedError(
            "EliminationVisitor::bvisit() not implemented for " + x.__str__()
        );
    }

    void EliminationVisitor::bvisit(const SymEngine::Symbol& x)
    {
        result_ = x.rcp_from_this();
    }

    void EliminationVisitor::bvisit(const SymEngine::Integer& x)
    {
        result_ = x.rcp_from_this();
    }

    void EliminationVisitor::bvisit(const SymEngine::Rational& x)
    {
        result_ = x.rcp_from_this();
    }

    void EliminationVisitor::bvisit(const SymEngine::Complex& x)
    {
        result_ = x.rcp_from_this();
    }

    void EliminationVisitor::bvisit(const SymEngine::Add& x)
    {
        SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
        // We check only `Basic` for each pair (`Basic` and `Number`) in the
        // dictionary of `Add`
        SymEngine::umap_basic_num d;
        for (const auto& p: x.get_dict()) {
            // Skip this pair if `Basic` was eliminated
            auto new_key = apply(p.first);
            if (!new_key.is_null()) SymEngine::Add::coef_dict_add_term(
                SymEngine::outArg(coef), d, p.second, new_key
            );
        }
        result_ = d.empty()
                ? SymEngine::RCP<const SymEngine::Basic>()
                : SymEngine::Add::from_dict(coef, std::move(d));
    }

    void EliminationVisitor::bvisit(const SymEngine::Mul& x)
    {
        SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
        // We check each pair (`Basic` and `Basic`) in the dictionary of `Mul`
        SymEngine::map_basic_basic d;
        for (const auto& p : x.get_dict()) {
            // Remove the whole `Mul` if the key will be eliminated
            auto new_key = apply(p.first);
            if (new_key.is_null()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
                return;
            }
            // The value in the pair (the exponent) is not allowed to be
            // eliminated completely
            auto new_value = apply(p.second);
            if (new_value.is_null()) throw SymEngine::SymEngineException(
                "EliminationVisitor::bvisit() does not allow to eliminate the exponent in a key-value pair of Mul."
            );
            SymEngine::Mul::dict_add_term_new(
                SymEngine::outArg(coef), d, new_value, new_key
            );
        }
        result_ = d.empty()
                ? SymEngine::RCP<const SymEngine::Basic>()
                : SymEngine::Mul::from_dict(coef, std::move(d));
    }

    void EliminationVisitor::bvisit(const SymEngine::Constant& x)
    {
        result_ = x.rcp_from_this();
    }

    void EliminationVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            eliminate_parameter<const NonElecFunction>(
                SymEngine::down_cast<const NonElecFunction&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            if (is_eliminable(op)) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                auto inner_state = op.get_inner_state();
                if (is_eliminable(*inner_state)) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    auto outer_state = op.get_outer_state();
                    if (is_eliminable(*outer_state)) {
                        result_ = SymEngine::RCP<const SymEngine::Basic>();
                    }
                    else {
                        result_ = x.rcp_from_this();
                    }
                }
            }
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            eliminate_one_arg_f<const CompositeFunction, const SymEngine::Basic>(
                op,
                op.get_inner(),
                [&](const SymEngine::RCP<const SymEngine::Basic>& inner)
                    -> SymEngine::RCP<const CompositeFunction>
                {
                    return SymEngine::make_rcp<const CompositeFunction>(
                        op.get_name(),
                        inner,
                        op.get_order()
                    );
                }
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            eliminate_one_arg_f<const ExchCorrEnergy, const SymEngine::Basic>(
                op,
                op.get_energy(),
                [&](const SymEngine::RCP<const SymEngine::Basic>& energy)
                    -> SymEngine::RCP<const ExchCorrEnergy>
                {
                    return SymEngine::make_rcp<const ExchCorrEnergy>(
                        op.get_name(),
                        op.get_state(),
                        op.get_overlap_distribution(),
                        op.get_weight(),
                        energy
                    );
                }
            );
        }
        else {
            throw SymEngine::NotImplementedError(
                "EliminationVisitor::bvisit() not implemented for FunctionSymbol " + x.__str__()
            );
        }
    }

    void EliminationVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        result_ = x.rcp_from_this();
    }

    void EliminationVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            eliminate_parameter<const OneElecDensity>(
                SymEngine::down_cast<const OneElecDensity&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            eliminate_parameter<const OneElecOperator>(
                SymEngine::down_cast<const OneElecOperator&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            if (is_eliminable(op)) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                auto state = op.get_state();
                if (is_eliminable(*state)) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    result_ = x.rcp_from_this();
                }
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            eliminate_one_arg_f<const ExchCorrPotential, const SymEngine::MatrixExpr>(
                op,
                op.get_potential(),
                [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& potential)
                    -> SymEngine::RCP<const ExchCorrPotential>
                {
                    return SymEngine::make_rcp<const ExchCorrPotential>(
                        op.get_name(),
                        op.get_state(),
                        op.get_overlap_distribution(),
                        op.get_weight(),
                        potential
                    );
                }
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            eliminate_one_arg_f<const TemporumOperator, const SymEngine::Basic>(
                op,
                op.get_target(),
                [&](const SymEngine::RCP<const SymEngine::Basic>& target)
                    -> SymEngine::RCP<const TemporumOperator>
                {
                    return SymEngine::make_rcp<const TemporumOperator>(
                        target, op.get_type()
                    );
                }
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            // We allow only the elimination of `x` as a whole
            eliminate_parameter<const TemporumOverlap>(
                SymEngine::down_cast<const TemporumOverlap&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const LagMultiplier>(x)) {
            eliminate_parameter<const LagMultiplier>(
                SymEngine::down_cast<const LagMultiplier&>(x)
            );
        }
        else {
            throw SymEngine::NotImplementedError(
                "EliminationVisitor::bvisit() not implemented for MatrixSymbol " + x.__str__()
            );
        }
    }

    void EliminationVisitor::bvisit(const SymEngine::Trace& x)
    {
        eliminate_one_arg_f<const SymEngine::Trace, const SymEngine::MatrixExpr>(
            x,
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(x.get_args()[0]),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Trace>
            {
                return SymEngine::make_rcp<const SymEngine::Trace>(arg);
            }
        );
    }

    void EliminationVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        eliminate_one_arg_f<const SymEngine::ConjugateMatrix, const SymEngine::MatrixExpr>(
            x,
            x.get_arg(),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::ConjugateMatrix>
            {
                return SymEngine::make_rcp<const SymEngine::ConjugateMatrix>(arg);
            }
        );
    }

    void EliminationVisitor::bvisit(const SymEngine::Transpose& x)
    {
        eliminate_one_arg_f<const SymEngine::Transpose, const SymEngine::MatrixExpr>(
            x,
            x.get_arg(),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Transpose>
            {
                return SymEngine::make_rcp<const SymEngine::Transpose>(arg);
            }
        );
    }

    void EliminationVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We check each argument of `MatrixAdd`
        SymEngine::vec_basic terms;
        for (auto arg: SymEngine::down_cast<const SymEngine::MatrixAdd&>(x).get_args()) {
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

    void EliminationVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // We check each argument of `MatrixMul`
        SymEngine::vec_basic factors;
        for (auto arg: SymEngine::down_cast<const SymEngine::MatrixMul&>(x).get_args()) {
            auto new_arg = apply(arg);
            if (new_arg.is_null()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
                return;
            }
            else {
                factors.push_back(new_arg);
            }
        }
        // Probably we do not need to check if factors is empty because the
        // above loop over arguments should be at least executed once
        if (factors.empty()) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            result_ = SymEngine::matrix_mul(factors);
        }
    }

    void EliminationVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        if (SymEngine::eq(*x.get_arg(), *parameter_)) {
            result_ = match_derivatives(x.get_symbols())
                    ? SymEngine::RCP<const SymEngine::Basic>()
                    : x.rcp_from_this();
        }
        else {
            result_ = x.rcp_from_this();
        }
    }
}
