#include <utility>

#include <symengine/number.h>
#include <symengine/symengine_exception.h>

#include "Tinned/operators/ConjugateTranspose.hpp"
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

#include "Tinned/visitors/EliminationVisitor.hpp"

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

    void EliminationVisitor::bvisit(const SymEngine::Number& x)
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
        // `SymEngine::Add::from_dict` will take care of empty `d`, that is
        // simply the coefficient
        result_ = SymEngine::Add::from_dict(coef, std::move(d));
    }

    void EliminationVisitor::bvisit(const SymEngine::Mul& x)
    {
        SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
        // We check each pair (`Basic` and `Basic`) in the dictionary of `Mul`
        SymEngine::map_basic_basic d;
        for (const auto& p: x.get_dict()) {
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
        // Probably we do not need to check if `d` is empty because the above
        // loop should either execute at least once or simply return
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
            result_ = x.rcp_from_this();
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            if (is_parameter_eliminable(op.get_inner_state())) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                if (is_parameter_eliminable(op.get_outer_state())) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    result_ = x.rcp_from_this();
                }
            }
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            eliminate_a_function(
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
            eliminate_a_function(
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
        if (SymEngine::is_a_sub<const PerturbedParameter>(x)) {
            eliminate_parameter(SymEngine::down_cast<const PerturbedParameter&>(x));
        }
        else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
            eliminate_a_function(
                op,
                std::bind(&construct_conjugate_transpose, std::placeholders::_1),
                op.get_arg()
            );
        }
        else if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            eliminate_parameter(SymEngine::down_cast<const OneElecDensity&>(x));
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            result_ = x.rcp_from_this();
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            if (is_parameter_eliminable(op.get_state())) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                result_ = x.rcp_from_this();
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            eliminate_a_function(
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
            eliminate_a_function(
                op,
                std::bind(&construct_dt_operator, std::placeholders::_1, op.get_type()),
                op.get_target()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            result_ = x.rcp_from_this();
        }
        else if (SymEngine::is_a_sub<const AdjointMap>(x)) {
            auto& op = SymEngine::down_cast<const AdjointMap&>(x);
            eliminate_a_function(
                op,
                std::bind(&construct_adjoint_map, std::placeholders::_1),
                op.get_x(),
                op.get_y()
            );
        }
        else if (SymEngine::is_a_sub<const ClusterConjHamiltonian>(x)) {
            auto& op = SymEngine::down_cast<const ClusterConjHamiltonian&>(x);
            eliminate_a_function(
                op,
                std::bind(&construct_cc_hamiltonian, std::placeholders::_1),
                op.get_cluster_operator(),
                op.get_hamiltonian()
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
        eliminate_a_function(
            x,
            std::bind(&construct_trace, std::placeholders::_1),
            x.get_args()[0]
        );
    }

    void EliminationVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        eliminate_a_function(
            x,
            std::bind(&construct_conjugate_matrix, std::placeholders::_1),
            x.get_arg()
        );
    }

    void EliminationVisitor::bvisit(const SymEngine::Transpose& x)
    {
        eliminate_a_function(
            x,
            std::bind(&construct_transpose, std::placeholders::_1),
            x.get_arg()
        );
    }

    void EliminationVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We check each argument of `MatrixAdd`
        SymEngine::vec_basic terms;
        for (auto arg: x.get_args()) {
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
        for (auto arg: x.get_args()) {
            auto new_arg = apply(arg);
            if (new_arg.is_null()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
                return;
            }
            else {
                factors.push_back(new_arg);
            }
        }
        // Probably we do not need to check if `factors` are empty because the
        // above loop should either execute at least once or simply return
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
