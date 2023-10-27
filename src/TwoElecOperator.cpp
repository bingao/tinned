#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/matrices/matrix_add.h>

#include "Tinned/TwoElecOperator.hpp"

namespace Tinned
{
    TwoElecOperator::TwoElecOperator(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const PertDependency& dependencies,
        const SymEngine::multiset_basic& derivative
    ) : SymEngine::MatrixSymbol(name),
        state_(state),
        dependencies_(dependencies),
        derivative_(derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t TwoElecOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *state_);
        for (auto& dep: dependencies_) {
            SymEngine::hash_combine(seed, *dep.first);
            SymEngine::hash_combine(seed, dep.second);
        }
        for (auto& p: derivative_) {
            SymEngine::hash_combine(seed, *p);
        }
        return seed;
    }

    bool TwoElecOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const TwoElecOperator>(o)) {
                auto& op = SymEngine::down_cast<const TwoElecOperator&>(o);
                return state_->__eq__(*op.state_)
                    && SymEngine::unified_eq(derivative_, op.derivative_)
                    && eq_dependency(dependencies_, op.dependencies_);
            }
        }
        return false;
    }

    int TwoElecOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TwoElecOperator>(o))
        int result = SymEngine::MatrixSymbol::compare(o);
        if (result == 0) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(o);
            result = state_->compare(*op.state_);
            if (result == 0) {
                result = SymEngine::unified_compare(derivative_, op.derivative_);
                return result == 0
                    ? SymEngine::ordered_compare(dependencies_, op.dependencies_)
                    : result;
            }
            else {
                return result;
            }
        }
        return result;
    }

    SymEngine::vec_basic TwoElecOperator::get_args() const
    {
        auto args = SymEngine::vec_basic({state_});
        auto deps = to_vec_basic(dependencies_);
        args.insert(args.end(), deps.begin(), deps.end());
        args.insert(args.end(), derivative_.begin(), derivative_.end());
        return args;
    }

    SymEngine::RCP<const SymEngine::Basic> TwoElecOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        // contr(g, D->diff(s))
        auto diff_state = SymEngine::rcp_dynamic_cast<const ElectronicState>(
            state_->diff(s)
        );
        auto op_diff_state = SymEngine::make_rcp<const TwoElecOperator>(
            get_name(),
            diff_state,
            dependencies_,
            derivative_
        );
        auto max_order = find_dependency(dependencies_, s);
        if (max_order > 0) {
            auto order = derivative_.count(s) + 1;
            if (order <= max_order) {
                // Return contr(g->diff(s), D) + contr(g, D->diff(s))
                auto derivative = derivative_;
                derivative.insert(s);
                return SymEngine::matrix_add({
                    SymEngine::make_rcp<const TwoElecOperator>(
                        get_name(),
                        state_,
                        dependencies_,
                        derivative
                    ),
                    op_diff_state
                });
            }
            else {
                return op_diff_state;
            }
        }
        else {
            return op_diff_state;
        }
    }
}
