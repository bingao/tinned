#include <symengine/symengine_exception.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

#include "Tinned/ReplaceVisitor.hpp"

namespace Tinned
{
    void ReplaceVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            // We only allow replacement of `NonElecFunction` as a whole, not
            // its arguments
            replace_whole<const NonElecFunction>(
                SymEngine::down_cast<const NonElecFunction&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            replace_one_arg_f<const CompositeFunction, const SymEngine::Basic>(
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
            replace_one_arg_f<const ExchCorrEnergy, const SymEngine::Basic>(
                op,
                op.get_energy(),
                [&](const SymEngine::RCP<const SymEngine::Basic>& energy)
                    -> SymEngine::RCP<const ExchCorrEnergy>
                {
                    return SymEngine::make_rcp<const ExchCorrEnergy>(op, energy);
                }
            );
        }
        else {
            SymEngine::MSubsVisitor::bvisit(x);
        }
    }

    void ReplaceVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        replace_whole<const SymEngine::ZeroMatrix>(x);
    }

    void ReplaceVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            replace_whole<const OneElecDensity>(
                SymEngine::down_cast<const OneElecDensity&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            replace_whole<const OneElecOperator>(
                SymEngine::down_cast<const OneElecOperator&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            replace_one_arg_f<const TwoElecOperator, const ElectronicState>(
                op,
                op.get_state(),
                [&](const SymEngine::RCP<const ElectronicState>& state)
                    -> SymEngine::RCP<const TwoElecOperator>
                {
                    return SymEngine::make_rcp<const TwoElecOperator>(
                        op.get_name(),
                        state,
                        op.get_dependencies(),
                        op.get_derivative()
                    );
                }
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            replace_one_arg_f<const ExchCorrPotential, const SymEngine::Basic>(
                op,
                op.get_potential(),
                [&](const SymEngine::RCP<const SymEngine::Basic>& potential)
                    -> SymEngine::RCP<const ExchCorrPotential>
                {
                    return SymEngine::make_rcp<const ExchCorrPotential>(op, potential);
                }
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            replace_one_arg_f<const TemporumOperator, const SymEngine::Basic>(
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
            replace_whole<const TemporumOverlap>(
                SymEngine::down_cast<const TemporumOverlap&>(x)
            );
        }
        else {
            SymEngine::MSubsVisitor::bvisit(x);
        }
    }

    void ReplaceVisitor::bvisit(const SymEngine::Trace& x)
    {
        replace_one_arg_f<const SymEngine::Trace, const SymEngine::MatrixExpr>(
            x,
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(x.get_args()[0]),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Trace>
            {
                return SymEngine::make_rcp<const SymEngine::Trace>(arg);
            }
        );
    }

    void ReplaceVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        replace_one_arg_f<const SymEngine::ConjugateMatrix, const SymEngine::MatrixExpr>(
            x,
            x.get_arg(),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::ConjugateMatrix>
            {
                return SymEngine::make_rcp<const SymEngine::ConjugateMatrix>(arg);
            }
        );
    }

    void ReplaceVisitor::bvisit(const SymEngine::Transpose& x)
    {
        replace_one_arg_f<const SymEngine::Transpose, const SymEngine::MatrixExpr>(
            x,
            x.get_arg(),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Transpose>
            {
                return SymEngine::make_rcp<const SymEngine::Transpose>(arg);
            }
        );
    }

    void ReplaceVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We first check each argument for replacement
        bool new_terms = false;
        SymEngine::vec_basic terms;
        for (auto arg: SymEngine::down_cast<const SymEngine::MatrixAdd&>(x).get_args()) {
            auto new_arg = apply(arg);
            if (SymEngine::neq(*arg, *new_arg)) new_terms = true;
            terms.push_back(new_arg);
        }
        SymEngine::RCP<const SymEngine::Basic> new_add;
        if (new_terms) {
            new_add = SymEngine::matrix_add(terms);
        }
        else {
            new_add = x.rcp_from_this();
        }
        // Next we check if the "new" `MatrixAdd` will be replaced as a whole
        replace_whole<const SymEngine::MatrixAdd>(
            SymEngine::down_cast<const SymEngine::MatrixAdd&>(*new_add)
        );
    }

    void ReplaceVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // We first check each argument for replacement
        bool new_factors = false;
        SymEngine::vec_basic factors;
        for (auto arg: SymEngine::down_cast<const SymEngine::MatrixMul&>(x).get_args()) {
            auto new_arg = apply(arg);
            if (SymEngine::neq(*arg, *new_arg)) new_factors = true;
            factors.push_back(new_arg);
        }
        SymEngine::RCP<const SymEngine::Basic> new_mul;
        if (new_factors) {
            new_mul = SymEngine::matrix_mul(factors);
        }
        else {
            new_mul = x.rcp_from_this();
        }
        // Next we check if the "new" `MatrixMul` will be replaced as a whole
        replace_whole<const SymEngine::MatrixMul>(
            SymEngine::down_cast<const SymEngine::MatrixMul&>(*new_mul)
        );
    }

    void ReplaceVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        // Because only `MatrixSymbol` can be used as the argument of
        // `MatrixDerivative`, we only need to check if `MatrixDerivative` will
        // be replaced as a whole
        replace_whole<const SymEngine::MatrixDerivative>(x);
    }
}
