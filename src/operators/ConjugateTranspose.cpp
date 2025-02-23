#include <cstddef>
#include <string>

#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/operators/ConjugateTranspose.hpp"
#include "Tinned/operators/ZeroOperator.hpp"
#include "Tinned/operators/TemporumOperator.hpp"

namespace Tinned
{
    ConjugateTranspose::ConjugateTranspose(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& arg
    ) : SymEngine::MatrixSymbol(std::string("Hermitian")),
        arg_(arg)
    {
        SYMENGINE_ASSIGN_TYPEID();
    }

    SymEngine::hash_t ConjugateTranspose::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *arg_);
        return seed;
    }

    bool ConjugateTranspose::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const ConjugateTranspose>(o)) {
            auto& op = SymEngine::down_cast<const ConjugateTranspose&>(o);
            return get_name()==op.get_name() && arg_->__eq__(*op.arg_);
        }
        return false;
    }

    int ConjugateTranspose::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ConjugateTranspose>(o))
        auto& op = SymEngine::down_cast<const ConjugateTranspose&>(o);
        if (get_name()==op.get_name()) {
            return arg_->compare(*op.arg_);
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic ConjugateTranspose::get_args() const
    {
        return SymEngine::vec_basic({arg_});
    }

    SymEngine::RCP<const SymEngine::Basic> ConjugateTranspose::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto diff_arg = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(
            arg_->diff(s)
        );
        if (diff_arg->__eq__(*make_zero_operator())) {
            return make_zero_operator();
        }
        else {
            return make_conjugate_transpose(diff_arg);
        }
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::Basic& x)
    {
        throw SymEngine::NotImplementedError(
            "ConjugateTransposeVisitor::bvisit() not implemented for " + x.__str__()
        );
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::MatrixExpr& x)
    {
        auto arg = SymEngine::rcp_static_cast<const SymEngine::MatrixExpr>(x.rcp_from_this());
        result_ = SymEngine::make_rcp<const ConjugateTranspose>(arg);
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::IdentityMatrix& x)
    {
        result_ = SymEngine::rcp_static_cast<const SymEngine::MatrixExpr>(x.rcp_from_this());
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        result_ = SymEngine::make_rcp<const SymEngine::ZeroMatrix>(x.ncols(), x.nrows());
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::DiagonalMatrix& x)
    {
        auto diag = x.get_container();
        SymEngine::vec_basic container(diag.size());
        for (std::size_t i=0; i<diag.size(); ++i)
            container[i] = SymEngine::conjugate(diag[i]);
        result_ = SymEngine::make_rcp<const SymEngine::DiagonalMatrix>(container);
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::ImmutableDenseMatrix& x)
    {
        SymEngine::vec_basic values(x.get_values().size());
        for (std::size_t i=0; i<x.nrows(); ++i)
            for (std::size_t j=0; j<x.ncols(); ++j)
                // [nrows, ncols] => [ncols, nrows]
                values[j*x.nrows()+i] = SymEngine::conjugate(x.get(i, j));
        result_ = SymEngine::make_rcp<const SymEngine::ImmutableDenseMatrix>(
            x.ncols(), x.nrows(), values
        );
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            TemporumType type = op.get_type()==TemporumType::Bra
                              ? TemporumType::Ket : TemporumType::Bra;
            result_ = make_dt_operator(apply(op.get_target()), type);
        }
        else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
            result_ = op.get_arg();
        }
        // For any Hermitian operator, its derivatives may not be Hermitian any more
        else {
            auto arg = SymEngine::rcp_static_cast<const SymEngine::MatrixExpr>(x.rcp_from_this());
            result_ = SymEngine::make_rcp<const ConjugateTranspose>(arg);
        }
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        result_ = SymEngine::transpose(x.get_arg());
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::Transpose& x)
    {
        result_ = SymEngine::conjugate_matrix(x.get_arg());
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        SymEngine::vec_basic terms;
        for (const auto& e: x.get_terms()) {
            e->accept(*this);
            terms.push_back(result_);
        }
        result_ = SymEngine::make_rcp<const SymEngine::MatrixAdd>(terms);
    }

    void ConjugateTransposeVisitor::bvisit(const SymEngine::HadamardProduct& x)
    {
        SymEngine::vec_basic factors;
        for (auto& e: x.get_factors()) {
            e->accept(*this);
            factors.push_back(result_);
        }
        result_ = SymEngine::make_rcp<const SymEngine::HadamardProduct>(factors);
    }
}

