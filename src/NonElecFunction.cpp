#include <symengine/integer.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/NonElecFunction.hpp"

namespace Tinned
{
    NonElecFunction::NonElecFunction(
        const std::string& name,
        const PertDependency& dependencies,
        const SymEngine::multiset_basic& derivatives
    ) : SymEngine::FunctionWrapper(name, SymEngine::vec_basic({})),
        dependencies_(dependencies),
        derivatives_(derivatives)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t NonElecFunction::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::FunctionWrapper::__hash__();
        hash_dependency(seed, dependencies_);
        for (auto& p: derivatives_) SymEngine::hash_combine(seed, *p);
        return seed;
    }

    bool NonElecFunction::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(o)) {
            auto& op = SymEngine::down_cast<const NonElecFunction&>(o);
            // We check the name, derivatives and perturbation dependencies
            return get_name() == op.get_name()
                && SymEngine::unified_eq(derivatives_, op.derivatives_)
                && eq_dependency(dependencies_, op.dependencies_);
        }
        return false;
    }

    int NonElecFunction::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const NonElecFunction>(o))
        auto& op = SymEngine::down_cast<const NonElecFunction&>(o);
        if (get_name() == op.get_name()) {
            int result = SymEngine::unified_compare(derivatives_, op.derivatives_);
            return result == 0
                ? SymEngine::ordered_compare(dependencies_, op.dependencies_)
                : result;
        }
        else {
            return get_name() < op.get_name() ? -1 : 1;
        }
    }

    //SymEngine::vec_basic NonElecFunction::get_args() const
    //{
    //    SymEngine::vec_basic args = to_vec_basic(dependencies_);
    //    args.insert(args.end(), derivatives_.begin(), derivatives_.end());
    //    return args;
    //}

    SymEngine::RCP<const SymEngine::Basic> NonElecFunction::create(
        const SymEngine::vec_basic &v
    ) const
    {
        throw SymEngine::NotImplementedError("NonElecFunction::create() not implemented");
    }

    SymEngine::RCP<const SymEngine::Number> NonElecFunction::eval(long bits) const
    {
        throw SymEngine::NotImplementedError("NonElecFunction::eval() not implemented");
    }

    SymEngine::RCP<const SymEngine::Basic> NonElecFunction::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto max_order = find_dependency(dependencies_, s);
        if (max_order > 0) {
            auto order = derivatives_.count(s) + 1;
            if (order <= max_order) {
                auto derivatives = derivatives_;
                derivatives.insert(s);
                return SymEngine::make_rcp<const NonElecFunction>(
                    get_name(),
                    dependencies_,
                    derivatives
                );
            }
            else {
                return SymEngine::integer(0);
            }
        }
        else {
            return SymEngine::integer(0);
        }
    }
}
