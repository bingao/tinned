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
        const SymEngine::multiset_basic& derivative
    ) : SymEngine::FunctionWrapper(name, SymEngine::vec_basic({})),
        dependencies_(dependencies),
        derivative_(derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t NonElecFunction::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::FunctionWrapper::__hash__();
        for (auto& dep: dependencies_) {
            SymEngine::hash_combine<const Perturbation>(seed, *dep.first);
            SymEngine::hash_combine<unsigned int>(seed, dep.second);
        }
        for (auto& p: derivative_) {
            SymEngine::hash_combine<SymEngine::Basic>(seed, *p);
        }
        return seed;
    }

    bool NonElecFunction::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::FunctionWrapper::__eq__(o)) {
            if (SymEngine::is_a_sub<const NonElecFunction>(o)) {
                const NonElecFunction& op = SymEngine::down_cast<const NonElecFunction &>(o);
                // First check the derivatives
                if (not SymEngine::unified_eq(derivative_, op.derivative_)) return false;
                // Secondly we check the perturbation dependencies
                return eq_dependency(dependencies_, op.dependencies_);
            }
        }
        return false;
    }

    int NonElecFunction::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const NonElecFunction>(o))
        int result = SymEngine::FunctionWrapper::compare(o);
        if (result == 0) {
            const NonElecFunction& op = SymEngine::down_cast<const NonElecFunction &>(o);
            result = SymEngine::unified_compare(derivative_, op.derivative_);
            if (result == 0) {
                return SymEngine::ordered_compare(dependencies_, op.dependencies_);
            }
            else {
                return result;
            }
        }
        return result;
    }

    SymEngine::vec_basic NonElecFunction::get_args() const
    {
        SymEngine::vec_basic args = to_vec_basic(dependencies_);
        args.insert(args.end(), derivative_.begin(), derivative_.end());
        return args;
    }

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
            auto order = derivative_.count(s) + 1;
            if (order <= max_order) {
                auto derivative = derivative_;
                derivative.insert(s);
                auto op = SymEngine::make_rcp<const NonElecFunction>(
                    SymEngine::FunctionWrapper::get_name(),
                    dependencies_,
                    derivative
                );
                return op;
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
