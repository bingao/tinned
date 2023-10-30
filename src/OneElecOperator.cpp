#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/OneElecOperator.hpp"

namespace Tinned
{
    OneElecOperator::OneElecOperator(
        const std::string& name,
        const PertDependency& dependencies,
        const SymEngine::multiset_basic& derivative
    ) : SymEngine::MatrixSymbol(name),
        dependencies_(dependencies),
        derivative_(derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t OneElecOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        for (auto& dep: dependencies_) {
            SymEngine::hash_combine(seed, *dep.first);
            SymEngine::hash_combine(seed, dep.second);
        }
        for (auto& p: derivative_) {
            SymEngine::hash_combine(seed, *p);
        }
        return seed;
    }

    bool OneElecOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const OneElecOperator>(o)) {
                auto& op = SymEngine::down_cast<const OneElecOperator&>(o);
                return SymEngine::unified_eq(derivative_, op.derivative_)
                    && eq_dependency(dependencies_, op.dependencies_);
            }
        }
        return false;
    }

    int OneElecOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const OneElecOperator>(o))
        int result = SymEngine::MatrixSymbol::compare(o);
        if (result == 0) {
            auto& op = SymEngine::down_cast<const OneElecOperator&>(o);
            result = SymEngine::unified_compare(derivative_, op.derivative_);
            return result == 0
                ? SymEngine::ordered_compare(dependencies_, op.dependencies_)
                : result;
        }
        return result;
    }

    //SymEngine::vec_basic OneElecOperator::get_args() const
    //{
    //    SymEngine::vec_basic args = to_vec_basic(dependencies_);
    //    args.insert(args.end(), derivative_.begin(), derivative_.end());
    //    return args;
    //}

    SymEngine::RCP<const SymEngine::Basic> OneElecOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto max_order = find_dependency(dependencies_, s);
        if (max_order > 0) {
            auto order = derivative_.count(s) + 1;
            if (order <= max_order) {
                auto derivative = derivative_;
                derivative.insert(s);
                return SymEngine::make_rcp<const OneElecOperator>(
                    get_name(),
                    dependencies_,
                    derivative
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
