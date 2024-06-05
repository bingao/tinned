#include <symengine/pow.h>

//#include "Tinned/Perturbation.hpp"
#include "Tinned/PerturbedParameter.hpp"
#include "Tinned/ConjugateTranspose.hpp"

#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecEnergy.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

#include "Tinned/AdjointMap.hpp"
#include "Tinned/ClusterConjHamiltonian.hpp"

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

    void FindAllVisitor::bvisit(const SymEngine::Number& x)
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
            auto outer = op.get_outer_state();
            // We first check electronic states
            if (find_only_name(*outer)) {
                find_only_name(*G->get_state());
            }
            // If electronic states are not those to find, we check
            // `TwoElecEnergy` and `TwoElecOperator`
            else {
                find_one_arg_f<const TwoElecEnergy, const TwoElecOperator>(
                    op,
                    [&](const TwoElecEnergy& op1, const TwoElecEnergy& op2) -> bool {
                        return op1.get_name()==op2.get_name()
                            && this->comp_2el_operator(
                                   *op1.get_2el_operator(), *op2.get_2el_operator()
                               )
                            && op1.get_outer_state()->get_name()
                               == op2.get_outer_state()->get_name();
                    },
                    G
                );
            }
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            find_one_arg_f<const CompositeFunction, const SymEngine::Basic>(
                op,
                [&](const CompositeFunction& op1, const CompositeFunction& op2) -> bool {
                    return op1.get_name()==op2.get_name()
                        && op1.get_inner()->__eq__(*op2.get_inner());
                },
                op.get_inner()
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            find_one_arg_f<const ExchCorrEnergy, const SymEngine::Basic>(
                op,
                [&](const ExchCorrEnergy& op1, const ExchCorrEnergy& op2) -> bool {
                    return op1.get_name()==op2.get_name()
                        && SymEngine::unified_eq(op1.get_args(), op2.get_args());
                },
                op.get_energy()
            );
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
        // We check only the name for the (perturbed) response parameter
        if (SymEngine::is_a_sub<const PerturbedParameter>(x)) {
            find_only_name(SymEngine::down_cast<const PerturbedParameter&>(x));
        }
        else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
            find_one_arg_f<const ConjugateTranspose, const SymEngine::MatrixExpr>(
                op,
                [&](const ConjugateTranspose& op1, const ConjugateTranspose& op2) -> bool
                {
                    FindAllVisitor v(op2.get_arg());
                    return !v.apply(op1.get_arg()).empty();
                },
                op.get_arg()
            );
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
            if (!find_with_condition<const TwoElecOperator>(
                op,
                [&](const TwoElecOperator& op1, const TwoElecOperator& op2) -> bool {
                    return this->comp_2el_operator(op1, op2);
                }
            )) find_only_name(*op.get_state());
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            find_one_arg_f<const ExchCorrPotential, const SymEngine::MatrixExpr>(
                op,
                [&](const ExchCorrPotential& op1, const ExchCorrPotential& op2) -> bool {
                    return op1.get_name()==op2.get_name()
                        && SymEngine::unified_eq(op1.get_args(), op2.get_args());
                },
                op.get_potential()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& dt = SymEngine::down_cast<const TemporumOperator&>(x);
            find_one_arg_f<const TemporumOperator, const SymEngine::MatrixExpr>(
                dt,
                // The strategy of comparing two `TemporumOperator` objects is
                // to make a `FindAllVisitor` with the target of the second
                // `TemporumOperator` object as the symbol to find, then apply
                // the visitor on the target of the first `TemporumOperator` object.
                [&](const TemporumOperator& op1, const TemporumOperator& op2) -> bool {
                    FindAllVisitor v(op2.get_target());
                    return !v.apply(op1.get_target()).empty();
                },
                dt.get_target()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            find_with_dependencies(SymEngine::down_cast<const TemporumOverlap&>(x));
        }
        else if (SymEngine::is_a_sub<const AdjointMap>(x)) {

        }
        else if (SymEngine::is_a_sub<const ClusterConjHamiltonian>(x)) {

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
        find_one_arg_f<const SymEngine::Trace, const SymEngine::Basic>(
            x,
            [&](const SymEngine::Trace& op1, const SymEngine::Trace& op2) -> bool {
                FindAllVisitor v(op2.get_args()[0]);
                return !v.apply(op1.get_args()[0]).empty();
            },
            x.get_args()[0]
        );
    }

    void FindAllVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        find_one_arg_f<const SymEngine::ConjugateMatrix, const SymEngine::MatrixExpr>(
            x,
            [&](const SymEngine::ConjugateMatrix& op1,
                const SymEngine::ConjugateMatrix& op2) -> bool {
                FindAllVisitor v(op2.get_arg());
                return !v.apply(op1.get_arg()).empty();
            },
            x.get_arg()
        );
    }

    void FindAllVisitor::bvisit(const SymEngine::Transpose& x)
    {
        find_one_arg_f<const SymEngine::Transpose, const SymEngine::MatrixExpr>(
            x,
            [&](const SymEngine::Transpose& op1, const SymEngine::Transpose& op2) -> bool
            {
                FindAllVisitor v(op2.get_arg());
                return !v.apply(op1.get_arg()).empty();
            },
            x.get_arg()
        );
    }

    void FindAllVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We first check if the whole `MatrixAdd` is that we are looking for,
        // if not, we check each argument of `MatrixAdd`
        if (!find_equivalence(x)) for (auto arg: x.get_args()) apply_(arg);
    }

    void FindAllVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // We first check if the whole `MatrixMul` is that we are looking for,
        // if not, we check each argument of `MatrixMul`
        if (!find_equivalence(x)) for (auto arg: x.get_args()) apply_(arg);
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
