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
        const SymEngine::multiset_basic& derivatives
    ) : SymEngine::MatrixSymbol(name),
        state_(state),
        dependencies_(dependencies),
        derivatives_(derivatives)
    {
        SYMENGINE_ASSERT(!is_zero_derivative(derivatives, dependencies))
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t TwoElecOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *state_);
        hash_dependency(seed, dependencies_);
        for (auto& p: derivatives_) SymEngine::hash_combine(seed, *p);
        return seed;
    }

    bool TwoElecOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const TwoElecOperator>(o)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(o);
            return get_name()==op.get_name()
                && state_->__eq__(*op.state_)
                && SymEngine::unified_eq(derivatives_, op.derivatives_)
                && eq_dependency(dependencies_, op.dependencies_);
        }
        return false;
    }

    int TwoElecOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TwoElecOperator>(o))
        auto& op = SymEngine::down_cast<const TwoElecOperator&>(o);
        if (get_name()==op.get_name()) {
            int result = state_->compare(*op.state_);
            if (result==0) {
                result = SymEngine::unified_compare(derivatives_, op.derivatives_);
                return result==0
                    ? SymEngine::ordered_compare(dependencies_, op.dependencies_)
                    : result;
            }
            else {
                return result;
            }
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic TwoElecOperator::get_args() const
    {
        return SymEngine::vec_basic({state_});
        //auto args = SymEngine::vec_basic({state_});
        //auto deps = dependency_to_vector(dependencies_);
        //args.insert(args.end(), deps.begin(), deps.end());
        //args.insert(args.end(), derivatives_.begin(), derivatives_.end());
        //return args;
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
            derivatives_
        );
        auto max_order = get_diff_order(s, dependencies_);
        if (max_order>0) {
            auto order = derivatives_.count(s) + 1;
            if (order<=max_order) {
                // Return contr(g->diff(s), D) + contr(g, D->diff(s))
                auto derivatives = derivatives_;
                derivatives.insert(s);
                return SymEngine::matrix_add({
                    SymEngine::make_rcp<const TwoElecOperator>(
                        get_name(),
                        state_,
                        dependencies_,
                        derivatives
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
