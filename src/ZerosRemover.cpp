#include <utility>

#include "Tinned/ZerosRemover.hpp"

namespace Tinned
{
    void ZerosRemover::bvisit(const SymEngine::Basic& x)
    {
        if (is_zero_quantity(x)) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            result_ = x.rcp_from_this();
        }
    }

    void ZerosRemover::bvisit(const SymEngine::Add& x)
    {
        SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
        // We check each pair (`Basic` and `Number`) in the dictionary of `Add`
        SymEngine::umap_basic_num d;
        for (const auto& p: x.get_dict()) {
            // Skip this pair if either `Basic` or `Number` is a zero quantity
            if (SymEngine::is_number_and_zero(*p.second)) continue;
            auto new_key = apply(p.first);
            if (!new_key.is_null()) SymEngine::Add::coef_dict_add_term(
                SymEngine::outArg(coef), d, p.second, new_key
            );
        }
        if (coef->is_zero() && d.empty()) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            result_ = SymEngine::Add::from_dict(coef, std::move(d));
        }
    }

    void ZerosRemover::bvisit(const SymEngine::Mul& x)
    {
        SymEngine::RCP<const SymEngine::Number> coef = x.get_coef();
        if (coef->is_zero()) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
            return;
        }
        // We check each pair (`Basic` and `Basic`) in the dictionary of `Mul`
        SymEngine::map_basic_basic d;
        for (const auto& p: x.get_dict()) {
            // Remove the whole `Mul` if the key is a zero quantity
            auto new_key = apply(p.first);
            if (new_key.is_null()) {
                result_ = SymEngine::RCP<const SymEngine::Basic>();
                return;
            }
            // Skip this pair if the value (the exponent) is a zero quantity
            auto new_value = apply(p.second);
            if (new_value.is_null()) continue;
            SymEngine::Mul::dict_add_term_new(
                SymEngine::outArg(coef), d, new_value, new_key
            );
        }
        // `SymEngine::Mul::from_dict` will take care of empty `d`
        result_ = SymEngine::Mul::from_dict(coef, std::move(d));
    }

    void ZerosRemover::bvisit(const SymEngine::Trace& x)
    {
        remove_one_arg_f<const SymEngine::Trace, const SymEngine::MatrixExpr>(
            x,
            SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(x.get_args()[0]),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::trace(arg);
            }
        );
    }

    void ZerosRemover::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        remove_one_arg_f<const SymEngine::ConjugateMatrix, const SymEngine::MatrixExpr>(
            x,
            x.get_arg(),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::conjugate_matrix(arg);
            }
        );
    }

    void ZerosRemover::bvisit(const SymEngine::Transpose& x)
    {
        remove_one_arg_f<const SymEngine::Transpose, const SymEngine::MatrixExpr>(
            x,
            x.get_arg(),
            [&](const SymEngine::RCP<const SymEngine::MatrixExpr>& arg)
                -> SymEngine::RCP<const SymEngine::Basic>
            {
                return SymEngine::transpose(arg);
            }
        );
    }

    void ZerosRemover::bvisit(const SymEngine::MatrixAdd& x)
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

    void ZerosRemover::bvisit(const SymEngine::MatrixMul& x)
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
        // Probably we do not need to check if `factors` are empty because the
        // above loop should either execute at least once or simply return
        if (factors.empty()) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            result_ = SymEngine::matrix_mul(factors);
        }
    }

    void ZerosRemover::bvisit(const SymEngine::MatrixDerivative& x)
    {
        auto new_arg = apply(x.get_arg());
        if (new_arg.is_null()) {
            result_ = SymEngine::RCP<const SymEngine::Basic>();
        }
        else {
            result_ = x.rcp_from_this();
        }
    }
}
