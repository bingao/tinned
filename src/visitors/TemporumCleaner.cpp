#include <cstddef>
#include <utility>

#include <symengine/number.h>
#include <symengine/symengine_assert.h>

#include "Tinned/operators/ConjugateTranspose.hpp"
#include "Tinned/operators/TemporumOperator.hpp"
#include "Tinned/operators/TemporumOverlap.hpp"

#include "Tinned/visitors/TemporumCleaner.hpp"

namespace Tinned
{
    void TemporumCleaner::bvisit(const SymEngine::Basic& x)
    {
        result_ = x.rcp_from_this();
    }

    void TemporumCleaner::bvisit(const SymEngine::Add& x)
    {
        SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
        // We check only `Basic` for each pair (`Basic` and `Number`) in the
        // dictionary of `Add`
        SymEngine::umap_basic_num d;
        for (const auto& p: x.get_dict()) {
            // Skip this pair if `Basic` becomes a zero quantity
            auto new_key = apply(p.first);
            if (!is_zero_quantity(new_key, threshold_))
                SymEngine::Add::coef_dict_add_term(
                    SymEngine::outArg(coef), d, p.second, new_key
                );
        }
        if (d.empty()) {
            result_ = SymEngine::zero;
        }
        else {
            result_ = SymEngine::Add::from_dict(coef, std::move(d));
        }
    }

    void TemporumCleaner::bvisit(const SymEngine::Mul& x)
    {
        SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
        // We check only the key for each pair (`Basic` and `Basic`) in the
        // dictionary of `Mul`
        SymEngine::map_basic_basic d;
        for (const auto& p: x.get_dict()) {
            // Set the whole `Mul` as zero if the key is zero quantity
            auto new_key = apply(p.first);
            if (is_zero_quantity(new_key, threshold_)) {
                result_ = SymEngine::zero;
                return;
            }
            // Exponent cannot be a `TemporumOperator` or `TemporumOverlap` object
            SymEngine::Mul::dict_add_term_new(
                SymEngine::outArg(coef), d, p.second, new_key
            );
        }
        if (d.empty()) {
            result_ = SymEngine::zero;
        }
        else {
            result_ = SymEngine::Mul::from_dict(coef, std::move(d));
        }
    }

    void TemporumCleaner::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
            clean_one_arg_f<const ConjugateTranspose, const SymEngine::MatrixExpr>(
                op,
                op.get_arg(),
                [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                    -> SymEngine::RCP<const SymEngine::Basic>
                {
                    return make_conjugate_transpose(arg);
                }
            );
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            // For unperturbed `TemporumOperator` objects, the function
            // `get_frequency()` will return zero frequency
            auto frequency = op.get_frequency();
            if (SymEngine::is_a_Number(*frequency)) {
                if (is_zero_number(
                    SymEngine::rcp_dynamic_cast<const SymEngine::Number>(frequency),
                    threshold_
                )) {
                    result_ = make_zero_operator();
                }
                else {
                    result_ = SymEngine::matrix_mul({frequency, op.get_target()});
                }
            }
            else {
                result_ = SymEngine::matrix_mul({frequency, op.get_target()});
            }
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOverlap&>(x);
            // `TemporumOverlap` will disappear if it is unperturbed or all
            // perturbations have zero frequencies
            for (std::size_t i=0; i<op.size(); ++i) {
                auto frequency = op.get_frequency(i);
                if (SymEngine::is_a_Number(*frequency)) {
                    if (!is_zero_number(
                        SymEngine::rcp_dynamic_cast<const SymEngine::Number>(frequency),
                        threshold_
                    )) {
                        result_ = x.rcp_from_this();
                        return;
                    }
                }
                else {
                    result_ = x.rcp_from_this();
                    return;
                }
            }
            result_ = make_zero_operator();
        }
        else {
            result_ = x.rcp_from_this();
        }
    }

    void TemporumCleaner::bvisit(const SymEngine::Trace& x)
    {
        clean_one_arg_f<const SymEngine::Trace, const SymEngine::MatrixExpr>(
            x,
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(x.get_args()[0]),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::trace(arg);
            }
        );
    }

    void TemporumCleaner::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        clean_one_arg_f<const SymEngine::ConjugateMatrix, const SymEngine::MatrixExpr>(
            x,
            x.get_arg(),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::conjugate_matrix(arg);
            }
        );
    }

    void TemporumCleaner::bvisit(const SymEngine::Transpose& x)
    {
        clean_one_arg_f<const SymEngine::Transpose, const SymEngine::MatrixExpr>(
            x,
            x.get_arg(),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::transpose(arg);
            }
        );
    }

    void TemporumCleaner::bvisit(const SymEngine::MatrixAdd& x)
    {
        // We check each argument of `MatrixAdd`
        SymEngine::vec_basic terms;
        for (auto arg: SymEngine::down_cast<const SymEngine::MatrixAdd&>(x).get_args()) {
            auto new_arg = apply(arg);
            if (!is_zero_quantity(new_arg, threshold_)) terms.push_back(new_arg);
        }
        if (terms.empty()) {
            result_ = make_zero_operator();
        }
        else {
            result_ = SymEngine::matrix_add(terms);
        }
    }

    void TemporumCleaner::bvisit(const SymEngine::MatrixMul& x)
    {
        // We check each argument of `MatrixMul`
        SymEngine::vec_basic factors;
        for (auto arg: SymEngine::down_cast<const SymEngine::MatrixMul&>(x).get_args()) {
            auto new_arg = apply(arg);
            if (is_zero_quantity(new_arg, threshold_)) {
                result_ = make_zero_operator();
                return;
            }
            else {
                factors.push_back(new_arg);
            }
        }
        // Probably we do not need to check if factors is empty because the
        // above loop over arguments should be at least executed once
        if (factors.empty()) {
            result_ = make_zero_operator();
        }
        else {
            result_ = SymEngine::matrix_mul(factors);
        }
    }

    void TemporumCleaner::bvisit(const SymEngine::MatrixDerivative& x)
    {
        // Although it is not forbidden, it is not a good idea to use a
        // `TemporumOperator` or `TemporumOverlap` object as the argument of
        // `SymEngine::MatrixDerivative`
        SymEngine::RCP<const SymEngine::Basic> arg = x.get_arg();
        if (SymEngine::is_a_sub<const TemporumOperator>(*arg) ||
            SymEngine::is_a_sub<const TemporumOverlap>(*arg)) {
            for (const auto& p: x.get_symbols()) {
                //FIXME: change stored `derivatives_` to the type
                //std::multiset<SymEngine::RCP<const SymEngine::Symbol> ?
                SYMENGINE_ASSERT(SymEngine::is_a_sub<const SymEngine::Symbol>(*p))
                auto s = SymEngine::rcp_dynamic_cast<const SymEngine::Symbol>(p);
                arg = arg->diff(s);
            }
            // Differentiated result is either an object of `TemporumOperator`
            // or `TemporumOverlap`, or zero
            result_ = apply(arg);
        }
        else {
            result_ = x.rcp_from_this();
        }
    }
}
