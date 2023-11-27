#include <symengine/add.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/TwoElecEnergy.hpp"

namespace Tinned
{
    TwoElecEnergy::TwoElecEnergy(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& inner,
        const SymEngine::RCP<const ElectronicState>& outer,
        const PertDependency& dependencies,
        const SymEngine::multiset_basic& derivative
    ) : SymEngine::FunctionWrapper(name, SymEngine::vec_basic({inner, outer})),
        dependencies_(dependencies),
        derivative_(derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t TwoElecEnergy::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::FunctionWrapper::__hash__();
        for (auto& dep: dependencies_) {
            SymEngine::hash_combine(seed, *dep.first);
            SymEngine::hash_combine(seed, dep.second);
        }
        for (auto& p: derivative_) {
            SymEngine::hash_combine(seed, *p);
        }
        return seed;
    }

    bool TwoElecEnergy::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const TwoElecEnergy>(o)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(o);
            return get_name() == op.get_name()
                && SymEngine::unified_eq(get_vec(), op.get_vec())
                && SymEngine::unified_eq(derivative_, op.derivative_)
                && eq_dependency(dependencies_, op.dependencies_);
        }
        return false;
    }

    int TwoElecEnergy::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TwoElecEnergy>(o))
        auto& op = SymEngine::down_cast<const TwoElecEnergy&>(o);
        if (get_name() == op.get_name()) {
            int result = SymEngine::unified_compare(get_vec(), op.get_vec());
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
        else {
            return get_name() < op.get_name() ? -1 : 1;
        }
    }

    //SymEngine::vec_basic TwoElecEnergy::get_args() const
    //{
    //    SymEngine::vec_basic args = to_vec_basic(dependencies_);
    //    args.insert(args.end(), derivative_.begin(), derivative_.end());
    //    return args;
    //}

    SymEngine::RCP<const SymEngine::Basic> TwoElecEnergy::create(
        const SymEngine::vec_basic &v
    ) const
    {
        throw SymEngine::NotImplementedError("TwoElecEnergy::create() not implemented");
    }

    SymEngine::RCP<const SymEngine::Number> TwoElecEnergy::eval(long bits) const
    {
        throw SymEngine::NotImplementedError("TwoElecEnergy::eval() not implemented");
    }

    SymEngine::RCP<const SymEngine::Basic> TwoElecEnergy::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        //FIXME: switch inner_ and outer_
        // tr(contr(g, inner_->diff(s)), outer_)
        auto diff_inner = SymEngine::rcp_dynamic_cast<const ElectronicState>(
            inner_->diff(s)
        );
        auto op_diff_inner = SymEngine::make_rcp<const TwoElecEnergy>(
            get_name(),
            diff_inner,
            outer_,
            dependencies_,
            derivative_
        );
        // tr(contr(g, inner_), outer_->diff(s))
        auto diff_outer = SymEngine::rcp_dynamic_cast<const ElectronicState>(
            outer_->diff(s)
        );
        auto op_diff_outer = SymEngine::make_rcp<const TwoElecEnergy>(
            get_name(),
            inner_,
            diff_outer,
            dependencies_,
            derivative_
        );
        auto op_diff_state = SymEngine::add(op_diff_inner, op_diff_outer);
        auto max_order = find_dependency(dependencies_, s);
        if (max_order > 0) {
            auto order = derivative_.count(s) + 1;
            if (order <= max_order) {
                auto derivative = derivative_;
                derivative.insert(s);
                return SymEngine::add(
                    // tr(contr(g->diff(s), inner_), outer_)
                    SymEngine::make_rcp<const TwoElecEnergy>(
                        get_name(),
                        inner_,
                        outer_,
                        dependencies_,
                        derivative
                    ),
                    op_diff_state
                );
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
