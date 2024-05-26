#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/OneElecOperator.hpp"
#include "Tinned/ZeroOperator.hpp"

namespace Tinned
{
    OneElecOperator::OneElecOperator(
        const std::string& name,
        const PertDependency& dependencies,
        const SymEngine::multiset_basic& derivatives
    ) : SymEngine::MatrixSymbol(name),
        dependencies_(dependencies),
        derivatives_(derivatives)
    {
        SYMENGINE_ASSERT(!is_zero_derivative(derivatives, dependencies))
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t OneElecOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        hash_dependency(seed, dependencies_);
        for (auto& p: derivatives_) SymEngine::hash_combine(seed, *p);
        return seed;
    }

    bool OneElecOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const OneElecOperator>(o)) {
            auto& op = SymEngine::down_cast<const OneElecOperator&>(o);
            return get_name()==op.get_name()
                && SymEngine::unified_eq(derivatives_, op.derivatives_)
                && eq_dependency(dependencies_, op.dependencies_);
        }
        return false;
    }

    int OneElecOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const OneElecOperator>(o))
        auto& op = SymEngine::down_cast<const OneElecOperator&>(o);
        if (get_name()==op.get_name()) {
            int result = SymEngine::unified_compare(derivatives_, op.derivatives_);
            return result==0
                ? SymEngine::ordered_compare(dependencies_, op.dependencies_)
                : result;
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    //SymEngine::vec_basic OneElecOperator::get_args() const
    //{
    //    SymEngine::vec_basic args = dependency_to_vector(dependencies_);
    //    args.insert(args.end(), derivatives_.begin(), derivatives_.end());
    //    return args;
    //}

    SymEngine::RCP<const SymEngine::Basic> OneElecOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto max_order = get_diff_order(s, dependencies_);
        if (max_order>0) {
            auto order = derivatives_.count(s) + 1;
            if (order<=max_order) {
                auto derivatives = derivatives_;
                derivatives.insert(s);
                return SymEngine::make_rcp<const OneElecOperator>(
                    get_name(),
                    dependencies_,
                    derivatives
                );
            }
            else {
                return make_zero_operator();
            }
        }
        else {
            return make_zero_operator();
        }
    }
}
