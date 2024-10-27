#include <symengine/pow.h>

#include "Tinned/PerturbedParameter.hpp"
#include "Tinned/ConjugateTranspose.hpp"

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

#include "Tinned/AdjointMap.hpp"
#include "Tinned/ClusterConjHamiltonian.hpp"

#include "Tinned/ExistAnyVisitor.hpp"

namespace Tinned
{
    void ExistAnyVisitor::bvisit(const SymEngine::Basic& x)
    {
        throw SymEngine::NotImplementedError(
            "ExistAnyVisitor::bvisit() not implemented for "+x.__str__()
        );
    }

    void ExistAnyVisitor::bvisit(const SymEngine::Symbol& x)
    {
        exist_any_equivalence(x);
    }

    void ExistAnyVisitor::bvisit(const SymEngine::Number& x)
    {
        exist_any_equivalence(x);
    }

    void ExistAnyVisitor::bvisit(const SymEngine::Add& x)
    {
        // We first check if the whole `Add` is that we are looking for
        if (!exist_any_equivalence(x)) {
            // We skip all pairs if the coefficient is what we are looking for
            if (!exist_any_equivalence(*x.get_coef())) {
                // Next we check each pair (`Basic` and `Number`) in the
                // dictionary of `Add`
                for (const auto& p: x.get_dict()) {
                    // Check if `Number` matches
                    if (exist_any_equivalence(*p.second)) break;
                    // Check if this pair is that we are looking for
                    if (exist_any_equivalence(*SymEngine::Add::from_dict(
                        SymEngine::zero, {{p.first, p.second}}
                    ))) break;
                    // Check if `Basic` matches
                    if (apply_(p.first)) break;
                }
            }
        }
    }

    void ExistAnyVisitor::bvisit(const SymEngine::Mul& x)
    {
        // We first check if the whole `Mul` is that we are looking for
        if (!exist_any_equivalence(x)) {
            // We skip all pairs if the coefficient is what we are looking for
            if (!exist_any_equivalence(*x.get_coef())) {
                // We check each pair (`Basic` and `Basic`) in the dictionary
                // of `Mul`
                for (const auto& p : x.get_dict()) {
                    // Check if this pair is that we are looking for
                    if (exist_any_equivalence(
                        *SymEngine::make_rcp<SymEngine::Pow>(p.first, p.second)
                    )) break;
                    // Check if the key and the value match
                    if (apply_(p.first)) break;
                    if (apply_(p.second)) break;
                }
            }
        }
    }

    void ExistAnyVisitor::bvisit(const SymEngine::Constant& x)
    {
        exist_any_equivalence(x);
    }

    void ExistAnyVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            exist_any_equivalence(x);
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
                if (!apply_(op.get_2el_operator())) apply_(op.get_outer_state());
            }
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
                apply_(op.get_inner());
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
                apply_(op.get_energy());
            }
        }
        else {
            throw SymEngine::NotImplementedError(
                "ExistAnyVisitor::bvisit() not implemented for FunctionSymbol "+x.__str__()
            );
        }
    }

    void ExistAnyVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        exist_any_equivalence(x);
    }

    void ExistAnyVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const PerturbedParameter>(x)) {
            exist_any_equivalence(x);
        }
        else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
                apply_(op.get_arg());
            }
        }
        else if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            exist_any_equivalence(x);
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            exist_any_equivalence(x);
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
                apply_(op.get_state());
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
                apply_(op.get_potential());
            }
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
                apply_(op.get_target());
            }
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            exist_any_equivalence(x);
        }
        else if (SymEngine::is_a_sub<const AdjointMap>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const AdjointMap&>(x);
                if (!apply_(op.get_y())) {
                    for (auto arg: op.get_x()) if (apply_(arg)) break;
                }
            }
        }
        else if (SymEngine::is_a_sub<const ClusterConjHamiltonian>(x)) {
            if (!exist_any_equivalence(x)) {
                auto& op = SymEngine::down_cast<const ClusterConjHamiltonian&>(x);
                if (!apply_(op.get_cluster_operator())) apply_(op.get_hamiltonian());
            }
        }
        else {
            throw SymEngine::NotImplementedError(
                "ExistAnyVisitor::bvisit() not implemented for MatrixSymbol "+x.__str__()
            );
        }
    }

    void ExistAnyVisitor::bvisit(const SymEngine::Trace& x)
    {
        if (!exist_any_equivalence(x)) apply_(x.get_args()[0]);
    }

    void ExistAnyVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        if (!exist_any_equivalence(x)) apply_(x.get_arg());
    }

    void ExistAnyVisitor::bvisit(const SymEngine::Transpose& x)
    {
        if (!exist_any_equivalence(x)) apply_(x.get_arg());
    }

    void ExistAnyVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        if (!exist_any_equivalence(x)) {
            for (auto arg: x.get_args()) if (apply_(arg)) break;
        }
    }

    void ExistAnyVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        if (!exist_any_equivalence(x)) {
            for (auto arg: x.get_args()) if (apply_(arg)) break;
        }
    }

    void ExistAnyVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        exist_any_equivalence(x);
    }
}
