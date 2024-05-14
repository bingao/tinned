#include "Tinned/Perturbation.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/LagMultiplier.hpp"

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

#include "Tinned/StateVector.hpp"
#include "Tinned/StateOperator.hpp"
#include "Tinned/AdjointMap.hpp"
#include "Tinned/ExpAdjointHamiltonian.hpp"

#include "Tinned/ZeroOperator.hpp"

#include "Tinned/ReplaceVisitor.hpp"

namespace Tinned
{
    void ReplaceVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            // We allow only replacement of `NonElecFunction` as a whole, not
            // its arguments
            replace_a_whole(SymEngine::down_cast<const NonElecFunction&>(x));
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            replace_a_function(
                op,
                [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
                {
                    return SymEngine::make_rcp<const TwoElecEnergy>(
                        SymEngine::rcp_dynamic_cast<const TwoElecOperator>(args[0]),
                        SymEngine::rcp_dynamic_cast<const ElectronicState>(args[1])
                    );
                },
                op.get_2el_operator(),
                op.get_outer_state()
            );
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            replace_a_function(
                op,
                [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
                {
                    return SymEngine::make_rcp<const CompositeFunction>(
                        op.get_name(), args[0], op.get_order()
                    );
                },
                op.get_inner()
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            // We also need to check the replacement of grid weight, state,
            // generalized overlap distribution
            replace_a_function(
                op,
                [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
                {
                    return SymEngine::make_rcp<const ExchCorrEnergy>(
                        op.get_name(),
                        SymEngine::rcp_dynamic_cast<const ElectronicState>(args[0]),
                        SymEngine::rcp_dynamic_cast<const OneElecOperator>(args[1]),
                        SymEngine::rcp_dynamic_cast<const NonElecFunction>(args[2]),
                        args[3]
                    );
                },
                op.get_state(),
                op.get_overlap_distribution(),
                op.get_weight(),
                op.get_energy()
            );
        }
        else {
            SymEngine::MSubsVisitor::bvisit(x);
        }
    }

    void ReplaceVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        replace_a_whole(x);
    }

    void ReplaceVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const LagMultiplier>(x)) {
            replace_a_whole(SymEngine::down_cast<const LagMultiplier&>(x));
        }
        else if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            replace_a_whole(SymEngine::down_cast<const OneElecDensity&>(x));
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            replace_a_whole(SymEngine::down_cast<const OneElecOperator&>(x));
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            replace_a_function(
                op,
                [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
                {
                    return SymEngine::make_rcp<const TwoElecOperator>(
                        op.get_name(),
                        SymEngine::rcp_dynamic_cast<const ElectronicState>(args[0]),
                        op.get_dependencies(),
                        op.get_derivatives()
                    );
                },
                op.get_state()
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            // We also need to check the replacement of grid weight, state,
            // generalized overlap distribution
            replace_a_function(
                op,
                [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
                {
                    return SymEngine::make_rcp<const ExchCorrPotential>(
                        op.get_name(),
                        SymEngine::rcp_dynamic_cast<const ElectronicState>(args[0]),
                        SymEngine::rcp_dynamic_cast<const OneElecOperator>(args[1]),
                        SymEngine::rcp_dynamic_cast<const NonElecFunction>(args[2]),
                        SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[3])
                    );
                },
                op.get_state(),
                op.get_overlap_distribution(),
                op.get_weight(),
                op.get_potential()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            replace_a_function(
                op,
                [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
                {
                    return SymEngine::make_rcp<const TemporumOperator>(
                        args[0], op.get_type()
                    );
                },
                op.get_target()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            replace_a_whole(SymEngine::down_cast<const TemporumOverlap&>(x));
        }
        else if (SymEngine::is_a_sub<const ZeroOperator>(x)) {
            replace_a_whole(SymEngine::down_cast<const ZeroOperator&>(x));
        }
        else {
            SymEngine::MSubsVisitor::bvisit(x);
        }
    }

    void ReplaceVisitor::bvisit(const SymEngine::Trace& x)
    {
        replace_a_function(
            x,
            [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::trace(
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0])
                );
            },
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(x.get_args()[0])
        );
    }

    void ReplaceVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        replace_a_function(
            x,
            [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::conjugate_matrix(
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0])
                );
            },
            x.get_arg()
        );
    }

    void ReplaceVisitor::bvisit(const SymEngine::Transpose& x)
    {
        replace_a_function(
            x,
            [&](const SymEngine::vec_basic& args) -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::transpose(
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(args[0])
                );
            },
            x.get_arg()
        );
    }

    void ReplaceVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // First we check if `x` will be replaced as a whole
        if (!replace_a_whole(x)) {
            // We next check each argument for replacement
            bool new_terms = false;
            SymEngine::vec_basic terms;
            for (auto& arg: x.get_args()) {
                auto new_arg = apply(arg);
                if (SymEngine::neq(*arg, *new_arg)) new_terms = true;
                terms.push_back(new_arg);
            }
            if (new_terms) {
                result_ = SymEngine::matrix_add(terms);
            }
            else {
                result_ = x.rcp_from_this();
            }
        }
    }

    void ReplaceVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // First we check if `x` will be replaced as a whole
        if (!replace_a_whole(x)) {
            // Next, we check each argument for replacement
            bool new_factors = false;
            SymEngine::vec_basic factors;
            for (auto& arg: x.get_args()) {
                auto new_arg = apply(arg);
                if (SymEngine::neq(*arg, *new_arg)) new_factors = true;
                factors.push_back(new_arg);
            }
            if (new_factors) {
                result_ = SymEngine::matrix_mul(factors);
            }
            else {
                result_ = x.rcp_from_this();
            }
        }
    }

    void ReplaceVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        // `MatrixDerivative` represents derivatives of a `MatrixSymbol`
        // object, so we need only check if `MatrixDerivative` will be replaced
        // as a whole instead of replacing the `MatrixSymbol` object
        replace_a_whole(x);
    }
}
