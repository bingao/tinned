#include <utility>

#include <symengine/constants.h>
#include <symengine/pow.h>
#include <symengine/symengine_exception.h>

#include "Tinned/RemoveVisitor.hpp"

namespace Tinned
{
    void RemoveIfVisitor::bvisit(const SymEngine::Basic& x)
    {
        throw SymEngine::NotImplementedError(
            "RemoveIfVisitor::bvisit() not implemented for " + x.__str__()
        );
    }

    void RemoveIfVisitor::bvisit(const SymEngine::Symbol& x)
    {
        remove_if_symbol_like<const SymEngine::Symbol>(x);
    }

    void RemoveIfVisitor::bvisit(const SymEngine::Integer& x)
    {
        remove_if_symbol_like<const SymEngine::Integer>(x);
    }

    void RemoveIfVisitor::bvisit(const SymEngine::Rational& x)
    {
        remove_if_symbol_like<const SymEngine::Rational>(x);
    }

    void RemoveIfVisitor::bvisit(const SymEngine::Complex& x)
    {
        remove_if_symbol_like<const SymEngine::Complex>(x);
    }

    void RemoveIfVisitor::bvisit(const SymEngine::Add& x)
    {
        // We first check if `Add` will be removed as a whole
        if (condition_(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            SymEngine::umap_basic_num d;
            // First we check if the coefficient will be removed
            SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
            if (condition_(*coef)) coef = SymEngine::zero;
            // Next we check each pair (`Basic` and `Number`) in the dictionary
            // of `Add`
            for (const auto& p: x.get_dict()) {
                // Skip if this pair will be removed as a whole
                if (condition_(*SymEngine::Add::from_dict(
                    SymEngine::zero, {{p.first, p.second}}
                ))) continue;
                // Skip if `Basic` was removed
                auto new_key = apply(p.first);
                if (new_key.is_null()) continue;
                // Skip if `Number` will be removed
                if (condition_(*p.second)) continue;
                SymEngine::Add::coef_dict_add_term(
                    SymEngine::outArg(coef), d, p.second, new_key
                );
            }
            result_ = SymEngine::Add::from_dict(coef, std::move(d));
        }
    }

    void RemoveIfVisitor::bvisit(const SymEngine::Mul& x)
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
                SymEngine::map_basic_basic d;
                // We check each pair (`Basic` and `Basic`) in the dictionary
                // of `Mul`
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

    void RemoveIfVisitor::bvisit(const SymEngine::Constant& x)
    {
        remove_if_symbol_like<const SymEngine::Constant>(x);
    }

    void RemoveIfVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        // We don't allow for the removal of derivative symbols, but only check
        // if the `NonElecFunction` (or its derivative) can be removed as a whole
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            remove_if_symbol_like<const SymEngine::NonElecFunction>(
                SymEngine::down_cast<const NonElecFunction&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {

        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveIfVisitor::bvisit() not implemented for FunctionSymbol " + x.__str__()
            );
        }
    }

    void RemoveIfVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        remove_if_symbol_like<const SymEngine::ZeroMatrix>(x);
    }

    void RemoveIfVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            remove_if_symbol_like<const OneElecDensity>(
                SymEngine::down_cast<const OneElecDensity&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            remove_if_symbol_like<const OneElecOperator>(
                SymEngine::down_cast<const OneElecOperator&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            remove_if_one_arg_f<const TwoElecOperator>(
                x,
                [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
                {
                    return op.get_state();
                },
                op.get_dependencies(),
                op.get_derivative()
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {

        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            remove_if_one_arg_f<const TemporumOperator>(
                x,
                [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
                {
                    return op.get_target();
                },
                op.get_type()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            remove_if_symbol_like<const TemporumOverlap>(
                SymEngine::down_cast<const TemporumOverlap&>(x)
            );
        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveIfVisitor::bvisit() not implemented for MatrixSymbol " + x.__str__()
            );
        }
    }

    void RemoveIfVisitor::bvisit(const SymEngine::Trace& x)
    {
        remove_if_one_arg_f<const SymEngine::Trace>(
            x,
            [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
            {
                return SymEngine::down_cast<const SymEngine::Trace&>(x).get_args()[0];
            }
        );
    }

    void RemoveIfVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        remove_ifnot_one_arg_f<const SymEngine::ConjugateMatrix>(
            x,
            [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
            {
                return SymEngine::down_cast<const SymEngine::ConjugateMatrix&>(x).get_arg();
            }
        );
    }

    void RemoveIfVisitor::bvisit(const SymEngine::Transpose& x)
    {
        remove_ifnot_one_arg_f<const SymEngine::Transpose>(
            x,
            [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
            {
                return SymEngine::down_cast<const SymEngine::Transpose&>(x).get_arg();
            }
        );
    }

    void RemoveIfVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We first check if `MatrixAdd` will be removed as a whole
        if (condition_(x)) {
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

    void RemoveIfVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // We first check if `MatrixMul` will be removed as a whole
        if (condition_(x)) {
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

    void RemoveIfVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        // Because only `MatrixSymbol` can be used as the argument of
        // `MatrixDerivative`, we only need to check if `MatrixDerivative` will
        // be removed as a whole
        remove_if_symbol_like<const SymEngine::MatrixDerivative>(x);
    }

    void RemoveIfNotVisitor::bvisit(const SymEngine::Mul& x)
    {
        if (condition_(x)) {
            SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
            // Indicator for keeping the whole `Mul`
            bool kept = !condition_(*coef);
            SymEngine::map_basic_basic d;
            // We check each pair (`Basic` and `Basic`) in the dictionary
            // of `Mul`
            for (const auto& p : x.get_dict()) {
                // First check the whole factor, i.e. key^value
                auto factor = SymEngine::make_rcp<SymEngine::Pow>(p.first, p.second);
                if (condition_(*factor)) {
                    // Skip this factor if the key will not be kept
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
                // The whole factor will be kept
                else {
                    kept = true;
                    SymEngine::Mul::dict_add_term_new(
                        SymEngine::outArg(coef), d, p.second, p.first
                    );
                }
            }
            result_ = SymEngine::Mul::from_dict(coef, std::move(d));
        }
        // `Mul` will be kept as a whole
        else {
            result_ = x.rcp_from_this();
        }
    }

    void RemoveIfNotVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        // We don't allow for the removal of derivative symbols, but only check
        // if the `NonElecFunction` (or its derivative) can be removed as a whole
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            remove_if_symbol_like<const SymEngine::NonElecFunction>(
                SymEngine::down_cast<const NonElecFunction&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {

        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveIfNotVisitor::bvisit() not implemented for FunctionSymbol " + x.__str__()
            );
        }
    }

    void RemoveIfNotVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            remove_if_symbol_like<const OneElecDensity>(
                SymEngine::down_cast<const OneElecDensity&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            remove_if_symbol_like<const OneElecOperator>(
                SymEngine::down_cast<const OneElecOperator&>(x)
            );
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            remove_ifnot_one_arg_f<const TwoElecOperator>(
                x,
                [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
                {
                    return op.get_state();
                },
                op.get_dependencies(),
                op.get_derivative()
            );
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {

        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            remove_ifnot_one_arg_f<const TemporumOperator>(
                x,
                [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
                {
                    return op.get_target();
                },
                op.get_type()
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            remove_if_symbol_like<const TemporumOverlap>(
                SymEngine::down_cast<const TemporumOverlap&>(x)
            );
        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveIfNotVisitor::bvisit() not implemented for MatrixSymbol " + x.__str__()
            );
        }
    }

    void RemoveIfNotVisitor::bvisit(const SymEngine::Trace& x)
    {
        remove_ifnot_one_arg_f<const SymEngine::Trace>(
            x,
            [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
            {
                return SymEngine::down_cast<const SymEngine::Trace&>(x).get_args()[0];
            }
        );
    }

    void RemoveIfNotVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        remove_ifnot_one_arg_f<const SymEngine::ConjugateMatrix>(
            x,
            [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
            {
                return SymEngine::down_cast<const SymEngine::ConjugateMatrix&>(x).get_arg();
            }
        );
    }

    void RemoveIfNotVisitor::bvisit(const SymEngine::Transpose& x)
    {
        remove_ifnot_one_arg_f<const SymEngine::Transpose>(
            x,
            [=]() -> SymEngine::RCP<const SymEngine::MatrixExpr>
            {
                return SymEngine::down_cast<const SymEngine::Transpose&>(x).get_arg();
            }
        );
    }

    void RemoveIfNotVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        // If `MatrixMul` will not be kept as whole, we then check if its
        // factors will be kept
        if (condition_(x)) {
            // `factors` will be used to construct the product that will be
            // removed. The first factor is -1 for substraction, see
            // information below.
            SymEngine::vec_basic factors = SymEngine::vec_basic({SymEngine::minus_one});
            // Indicates if there is factor(s) kept
            bool factors_kept = false;
            for (auto arg: SymEngine::down_cast<const SymEngine::MatrixMul&>(x).get_args()) {
                auto new_arg = apply(arg);
                if (new_arg.is_null()) {
                    // This factor does not match any given symbols, but we
                    // save it in case that there will be factor(s) kept
                    factors.push_back(arg);
                }
                else {
                    // `MatrixMul` will be kept as a whole if one of its
                    // factors match any given symbols and is kept
                    if (SymEngine::eq(*arg, *new_arg)) {
                        result_ = x.rcp_from_this();
                        return;
                    }
                    else {
                        // Suppose `MatrixMul` is A*B*C*... = (Ak+Ar)*B*C*...,
                        // where Ak will be kept and Ar will be removed. The
                        // result after removal will be Ak*B*C*..., and we save
                        // Ar = A-Ak. The result can also be computed as
                        // A*B*C*... - Ar*B*C*...
                        factors.push_back(SymEngine::matrix_add(
                            SymEngine::vec_basic({
                                arg,
                                SymEngine::matrix_mul({SymEngine::minus_one, new_arg})
                            })
                        ));
                        factors_kept = true;
                    }
                }
            }
            // As aforementioned, when there are factors partially kept, the
            // result can be computed as A*B*C*...*R*S*T*... - Ar*Br*Cr*...*R*S*T*...,
            // where Ar, Br, Cr, ... are parts that are removed
            if (factors_kept) {
                result_ = SymEngine::matrix_add(
                    SymEngine::vec_basic({
                        x.rcp_from_this(),
                        SymEngine::matrix_mul(factors)
                    })
                );
            }
            // `MatrixMul` will be removed if all its factors are null after
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
