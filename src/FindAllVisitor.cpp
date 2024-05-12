#include <iostream>

#include <symengine/pow.h>

//#include "Tinned/Perturbation.hpp"
#include "Tinned/LagMultiplier.hpp"

#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecEnergy.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

#include "Tinned/StateVector.hpp"
#include "Tinned/StateOperator.hpp"
#include "Tinned/AdjointMap.hpp"
#include "Tinned/ExpAdjointHamiltonian.hpp"

#include "Tinned/ZeroOperator.hpp"

#include "Tinned/FindAllVisitor.hpp"

namespace Tinned
{
    void FindAllVisitor::bvisit(const SymEngine::Basic& x)
    {
        throw SymEngine::NotImplementedError(
            "FindAllVisitor::bvisit() not implemented for "+x.__str__()
        );
    }

    void FindAllVisitor::bvisit(const SymEngine::Symbol& x)
    {
        find_equivalence(x);
    }

    void FindAllVisitor::bvisit(const SymEngine::Integer& x)
    {
        find_equivalence(x);
    }

    void FindAllVisitor::bvisit(const SymEngine::Rational& x)
    {
        find_equivalence(x);
    }

    void FindAllVisitor::bvisit(const SymEngine::Complex& x)
    {
        find_equivalence(x);
    }

    void FindAllVisitor::bvisit(const SymEngine::Add& x)
    {
        // We first check if the whole `Add` is that we are looking for
        if (!find_equivalence(x)) {
            // We skip all pairs if the coefficient is what we are looking for
            if (!find_equivalence(*x.get_coef())) {
                // Next we check each pair (`Basic` and `Number`) in the
                // dictionary of `Add`
                for (const auto& p: x.get_dict()) {
                    // Check if `Number` matches
                    if (find_equivalence(*p.second)) continue;
                    // Check if this pair is that we are looking for
                    if (find_equivalence(*SymEngine::Add::from_dict(
                        SymEngine::zero, {{p.first, p.second}}
                    ))) continue;
                    // Check if `Basic` matches
                    apply_(p.first);
                }
            }
        }
    }

    void FindAllVisitor::bvisit(const SymEngine::Mul& x)
    {
        // We first check if the whole `Mul` is that we are looking for
        if (!find_equivalence(x)) {
            // We skip all pairs if the coefficient is what we are looking for
            if (!find_equivalence(*x.get_coef())) {
                // We check each pair (`Basic` and `Basic`) in the dictionary
                // of `Mul`
                for (const auto& p : x.get_dict()) {
                    // Check if this pair is that we are looking for
                    if (find_equivalence(
                        *SymEngine::make_rcp<SymEngine::Pow>(p.first, p.second)
                    )) continue;
                    // Check if the key and the value match
                    apply_(p.first);
                    apply_(p.second);
                }
            }
        }
    }

    void FindAllVisitor::bvisit(const SymEngine::Constant& x)
    {
        find_equivalence(x);
    }

    void FindAllVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            find_with_dependencies(SymEngine::down_cast<const NonElecFunction&>(x));
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            auto G = op.get_2el_operator();
            if (!find_2el_operator(*G)) {
                auto inner = G->get_state();
                auto outer = op.get_outer_state();
                // We need to find all different states
                if (find_only_name(*inner)) find_only_name(*outer);
                // For `TwoElecEnergy`, we check its name, its two-electron
                // operator's name and dependencies, and its inner and outer
                // density matrices' names
                else if (SymEngine::is_a_sub<const TwoElecEnergy>(*symbol_)) {
                    auto s = SymEngine::rcp_dynamic_cast<const TwoElecEnergy>(symbol_);
                    auto s_G = s->get_2el_operator();
                    if (op.get_name()==s->get_name() &&
                        G->get_name()==s_G->get_name() &&
                        eq_dependency(G->get_dependencies(), s_G->get_dependencies()) &&
                        inner->get_name()==s_G->get_state()->get_name() &&
                        outer->get_name()==s->get_outer_state()->get_name())
                        result_.insert(x.rcp_from_this());
                }
            }
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            auto inner = op.get_inner();
            // We first check if the composite function is that we are finding,
            // by comparing its name and inner function
            if (SymEngine::is_a_sub<const CompositeFunction>(*symbol_)) {
                auto s = SymEngine::rcp_dynamic_cast<const CompositeFunction>(symbol_);
                if (op.get_name()==s->get_name() && inner->__eq__(*s->get_inner()))
                    result_.insert(x.rcp_from_this());
            }
            // If the composite function is not we are looking for, we check
            // its inner function
            else {
                apply_(inner);
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            // We first check if `ExchCorrEnergy` is that we are finding, by
            // comparing its name and arguments
            if (SymEngine::is_a_sub<const ExchCorrEnergy>(*symbol_)) {
                auto s = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(symbol_);
                if (op.get_name()==s->get_name() &&
                    SymEngine::unified_eq(op.get_args(), s->get_args())) {
                    result_.insert(x.rcp_from_this());
                }
            }
            // If the `ExchCorrEnergy` is not we are looking for, we check
            // its XC energy or derivatives
            else {
                apply_(op.get_energy());
            }
        }
        else {
            throw SymEngine::NotImplementedError(
                "FindAllVisitor::bvisit() not implemented for FunctionSymbol "+x.__str__()
            );
        }
    }

    void FindAllVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        find_equivalence(x);
    }

    void FindAllVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        // We check only the name for the Lagrangian multiplier
        if (SymEngine::is_a_sub<const LagMultiplier>(x)) {
            find_only_name(SymEngine::down_cast<const LagMultiplier&>(x));
        }
        // We check only the name for one-electron spin-orbital density matrix
        else if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            find_only_name(SymEngine::down_cast<const OneElecDensity&>(x));
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            find_with_dependencies(SymEngine::down_cast<const OneElecOperator&>(x));
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            if (!find_2el_operator(op)) find_only_name(*op.get_state());
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            // We first check if `ExchCorrPotential` is that we are finding, by
            // comparing its name and arguments
            if (SymEngine::is_a_sub<const ExchCorrPotential>(*symbol_)) {
                auto s = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(symbol_);
                if (op.get_name()==s->get_name() &&
                    SymEngine::unified_eq(op.get_args(), s->get_args())) {
                    result_.insert(x.rcp_from_this());
                }
            }
            // If the `ExchCorrPotential` is not we are looking for, we check
            // its XC potential operator or derivatives
            else {
                apply_(op.get_potential());
            }
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& dt = SymEngine::down_cast<const TemporumOperator&>(x);
            find_one_arg_f<const TemporumOperator, const SymEngine::Basic>(
                dt,
                dt.get_target(),
                [&](const SymEngine::RCP<const TemporumOperator>& op) -> FindAllVisitor
                {
                    return FindAllVisitor(op->get_target());
                }
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            find_with_dependencies(SymEngine::down_cast<const TemporumOverlap&>(x));
        }
        else if (SymEngine::is_a_sub<const ZeroOperator>(x)) {
            find_equivalence(SymEngine::down_cast<const ZeroOperator&>(x));
        }
        else {
            throw SymEngine::NotImplementedError(
                "FindAllVisitor::bvisit() not implemented for MatrixSymbol "+x.__str__()
            );
        }
    }

    // For `Trace`, `ConjugateMatrix` and `Transpose`, we follow the same
    // procedure as `TemporumOperator`
    void FindAllVisitor::bvisit(const SymEngine::Trace& x)
    {
        find_one_arg_f<const SymEngine::Trace, const SymEngine::MatrixExpr>(
            x,
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(x.get_args()[0]),
            [&](const SymEngine::RCP<const SymEngine::Trace>& op) -> FindAllVisitor
            {
                return FindAllVisitor(op->get_args()[0]);
            }
        );
    }

    void FindAllVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        find_one_arg_f<const SymEngine::ConjugateMatrix, const SymEngine::MatrixExpr>(
            x,
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(x.get_arg()),
            [&](const SymEngine::RCP<const SymEngine::ConjugateMatrix>& op) -> FindAllVisitor
            {
                return FindAllVisitor(op->get_arg());
            }
        );
    }

    void FindAllVisitor::bvisit(const SymEngine::Transpose& x)
    {
        find_one_arg_f<const SymEngine::Transpose, const SymEngine::MatrixExpr>(
            x,
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(x.get_arg()),
            [&](const SymEngine::RCP<const SymEngine::Transpose>& op) -> FindAllVisitor
            {
                return FindAllVisitor(op->get_arg());
            }
        );
    }

    void FindAllVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We first check if the whole `MatrixAdd` is that we are looking for
        if (!find_equivalence(x)) {
            // Next we check each argument of `MatrixAdd`
            for (auto arg: SymEngine::down_cast<const SymEngine::MatrixAdd&>(x).get_args())
                apply_(arg);
        }
    }

    void FindAllVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // We first check if the whole `MatrixMul` is that we are looking for
        if (!find_equivalence(x)) {
            // Next we check each argument of `MatrixMul`
            for (auto arg: SymEngine::down_cast<const SymEngine::MatrixMul&>(x).get_args())
                apply_(arg);
        }
    }

    void FindAllVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        // `MatrixDerivative` represents derivatives of a `MatrixSymbol`
        // object, so we need only check if its this object is that we are
        // looking for
        auto arg = x.get_arg();
        if (arg->__eq__(*symbol_)) result_.insert(x.rcp_from_this());
    }
}
