#include <cstddef>
#include <algorithm>
#include <string>

#include <symengine/matrices/matrix_add.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/operators/AdjointMap.hpp"
#include "Tinned/operators/ZeroOperator.hpp"

namespace Tinned
{
    AdjointMap::AdjointMap(
        const SymEngine::vec_basic& x,
        const SymEngine::RCP<const SymEngine::Basic>& y
    ) : SymEngine::MatrixSymbol(std::string("ad"))
    {
        x_ = x;
        if (x_.size()>1) std::sort(x_.begin(), x_.end(), SymEngine::RCPBasicKeyLess());
        y_ = y;
        SYMENGINE_ASSIGN_TYPEID()
    }

    AdjointMap::AdjointMap(
        const AdjointMap& other,
        const SymEngine::RCP<const SymEngine::Basic>& x
    ) : SymEngine::MatrixSymbol(other.get_name())
    {
        x_ = other.x_;
        x_.push_back(x);
        if (x_.size()>1) std::sort(x_.begin(), x_.end(), SymEngine::RCPBasicKeyLess());
        y_ = other.y_;
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t AdjointMap::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        for (auto& op: x_) SymEngine::hash_combine(seed, *op);
        SymEngine::hash_combine(seed, *y_);
        return seed;
    }

    bool AdjointMap::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const AdjointMap>(o)) {
            auto& op = SymEngine::down_cast<const AdjointMap&>(o);
            return get_name()==op.get_name()
                && SymEngine::unified_eq(x_, op.x_)
                && y_->__eq__(*op.y_);
        }
        return false;
    }

    int AdjointMap::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const AdjointMap>(o))
        auto& op = SymEngine::down_cast<const AdjointMap&>(o);
        if (get_name()==op.get_name()) {
            int result = SymEngine::unified_compare(x_, op.x_);
            return result==0 ? y_->compare(*op.y_) : result;
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic AdjointMap::get_args() const
    {
        SymEngine::vec_basic args = x_;
        args.push_back(y_);
        return args;
    }

    SymEngine::RCP<const SymEngine::Basic> AdjointMap::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        SymEngine::vec_basic terms;
        // Differentiate each X
        for (std::size_t i=0; i<x_.size(); ++i) {
            auto diff_x = x_;
            diff_x[i] = x_[i]->diff(s);
            // Skip zero derivative
            if (diff_x[i]->__eq__(*make_zero_operator())) continue;
            // Sort X's
            if (diff_x.size()>1)
                std::sort(diff_x.begin(), diff_x.end(), SymEngine::RCPBasicKeyLess());
            terms.push_back(SymEngine::make_rcp<const AdjointMap>(diff_x, y_));
        }
        // Differentiate Y
        auto diff_y = y_->diff(s);
        if (!diff_y->__eq__(*make_zero_operator())) {
            terms.push_back(SymEngine::make_rcp<const AdjointMap>(x_, diff_y));
        }
        if (terms.empty()) {
            return make_zero_operator();
        }
        else {
            if (terms.size()==1) {
                return terms[0];
            }
            else {
                return SymEngine::matrix_add(terms);
            }
        }
    }
}
