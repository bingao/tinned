#include <utility>

#include <symengine/pow.h>
#include <symengine/symengine_exception.h>

#include "Tinned/RemoveVisitor.hpp"

namespace Tinned
{
    RemoveVisitor::RemoveVisitor(const SymEngine::vec_basic& symbols, bool equivalence)
        : symbols_(symbols)
    {
        if (equivalence) {
            to_remove_ = [=](const SymEngine::Basic& x) -> bool
            {
                return this->eq_check(x);
            };
        }
        else {
            to_remove_ = [=](const SymEngine::Basic& x) -> bool
            {
                return this->neq_check(x);
            };
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Basic& x)
    {
        throw SymEngine::NotImplementedError(
            "RemoveVisitor::bvisit() not implemented for " + x.__str__()
        );
    }

    // Each upper level check if a zero matrix/symbol returned from lower level

    void RemoveVisitor::bvisit(const SymEngine::Symbol& x)
    {
        result_ = to_remove_(x)
            ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
    }

    void RemoveVisitor::bvisit(const SymEngine::Integer& x)
    {
        result_ = to_remove_(x)
            ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
    }

    void RemoveVisitor::bvisit(const SymEngine::Rational& x)
    {
        result_ = to_remove_(x)
            ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
    }

    void RemoveVisitor::bvisit(const SymEngine::Complex& x)
    {
        result_ = to_remove_(x)
            ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
    }

    void RemoveVisitor::bvisit(const SymEngine::Add& x)
    {
        // We first check if `Add` will be removed as a whole
        if (to_remove_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            SymEngine::umap_basic_num d;
            // First we check if the coefficient will be removed
            SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
            if (to_remove_(*coef)) coef = SymEngine::zero;
            // Next we check each pair (`Basic` and `Number`) in the dictionary
            // of `Add`
            for (const auto& p: x.get_dict()) {
                // Skip if this pair will be removed as a whole
                if (to_remove_(*SymEngine::Add::from_dict(
                    SymEngine::zero, {{p.first, p.second}}
                ))) continue;
                // Skip if `Basic` was removed
                auto new_key = apply(p.first);
                if (new_key.is_null()) continue;
                // Skip if `Number` will be removed
                if (to_remove_(*p.second)) continue;
                SymEngine::Add::coef_dict_add_term(
                    SymEngine::outArg(coef), d, p.second, new_key
                );
            }
            result_ = SymEngine::Add::from_dict(coef, std::move(d));
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Mul& x)
    {
        // We first check if `Mul` will be removed as a whole
        if (to_remove_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            // We remove the whole `Mul` if its coefficient will be removed
            SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
            if (to_remove_(*coef)) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                SymEngine::map_basic_basic d;
                // We check each pair (`Basic` and `Basic`) in the dictionary
                // of `Mul`
                for (const auto& p : x.get_dict()) {
                    // We remove the whole `Mul` if the pair will be removed
                    auto factor = SymEngine::make_rcp<SymEngine::Pow>(p.first, p.second);
                    if (to_remove_(*factor)) {
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
                    if (new_key.is_null()) {
                        throw SymEngine::SymEngineException(
                            "Removing the exponent in a key-value pair of Mul is not allowed."
                        );
                    }
                    else {
                        SymEngine::Mul::dict_add_term_new(
                            SymEngine::outArg(coef), d, new_value, new_key
                        );
                    }
                }
                result_ = SymEngine::Mul::from_dict(coef, std::move(d));
            }
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Constant& x)
    {
        result_ = to_remove_(x)
            ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
    }

    void RemoveVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        // We don't allow for the removal of derivative symbols, but only check
        // if the `NonElecFunction` (or its derivative) can be removed as a whole
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            result_ = to_remove_(x)
                ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {

        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveVisitor::bvisit() not implemented for FunctionSymbol " + x.__str__()
            );
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        result_ = to_remove_(x)
            ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            result_ = to_remove_(x)
                ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            result_ = to_remove_(x)
                ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            // We first check if `TwoElecOperator` will be removed
            if (to_remove_(x)) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            // Next we check if its state will be removed
            else {
                auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
                auto new_state =  apply(op.get_state());
                if (new_state.is_null()) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    result_ = SymEngine::make_rcp<const TwoElecOperator>(
                        op.get_name(),
                        SymEngine::rcp_dynamic_cast<const ElectronicState>(new_state),
                        op.get_dependencies(),
                        op.get_derivative()
                    );
                }
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {

        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            // We first check if `TemporumOperator` will be removed
            if (to_remove_(x)) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            // Next we check if its target will be removed
            else {
                auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
                auto new_target = apply(op.get_target());
                if (new_target.is_null()) {
                    result_ = SymEngine::RCP<const SymEngine::Basic>();
                }
                else {
                    result_ = SymEngine::make_rcp<const TemporumOperator>(
                        new_target, op.get_type()
                    );
                }
            }
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            result_ = to_remove_(x)
                ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveVisitor::bvisit() not implemented for MatrixSymbol " + x.__str__()
            );
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Trace& x)
    {
        // We first check if `Trace` will be removed
        if (to_remove_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        // Next we check if its argument will be removed
        else {
            auto new_arg = apply(
                SymEngine::down_cast<const SymEngine::Trace&>(x).get_args()[0]
            );
            if (new_arg.is_null()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                result_ = SymEngine::make_rcp<const SymEngine::Trace>(
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(new_arg)
                );
            }
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        // We first check if `ConjugateMatrix` will be removed
        if (to_remove_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        // Next we check if its argument will be removed
        else {
            auto new_arg = apply(
                SymEngine::down_cast<const SymEngine::ConjugateMatrix&>(x).get_arg()
            );
            if (new_arg.is_null()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                result_ = SymEngine::make_rcp<const SymEngine::ConjugateMatrix>(
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(new_arg)
                );
            }
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Transpose& x)
    {
        // We first check if `Transpose` will be removed
        if (to_remove_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        // Next we check if its argument will be removed
        else {
            auto new_arg = apply(
                SymEngine::down_cast<const SymEngine::Transpose&>(x).get_arg()
            );
            if (new_arg.is_null()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
            }
            else {
                result_ = SymEngine::make_rcp<const SymEngine::Transpose>(
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(new_arg)
                );
            }
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We first check if `MatrixAdd` will be removed as a whole
        if (to_remove_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            // Next we check each argument of `MatrixAdd`
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
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // We first check if `MatrixMul` will be removed as a whole
        if (to_remove_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            // Next we check each argument of `MatrixMul`
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
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        // Because only `MatrixSymbol` can be used as the argument of
        // `MatrixDerivative`, we only need to check if `MatrixDerivative` will
        // be removed as a whole
        result_ = to_remove_(x)
            ? SymEngine::RCP<const SymEngine::Basic>() : x.rcp_from_this();
    }
}
