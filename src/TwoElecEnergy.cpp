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
        const SymEngine::multiset_basic& derivatives
    ) : SymEngine::FunctionWrapper(name, SymEngine::vec_basic({inner, outer}))
    {
        // For numerical evaulation, one may put fewer density matrices in the
        // innermost loop. But we do not know the number of components of each
        // perturbation in Tinned, we also do not know the number of components
        // of high order differentiation with respect to a same perturbation.
        // So we simply rearrange inner_ and outer_ according to their
        // `compare` function.
        if (inner->compare(*outer)<=0) {
            inner_ = inner;
            outer_ = outer;
        }
        else {
            inner_ = outer;
            outer_ = inner;
        }
        dependencies_ = dependencies;
        derivatives_ = derivatives;
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t TwoElecEnergy::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::FunctionWrapper::__hash__();
        hash_dependency(seed, dependencies_);
        for (auto& p: derivatives_) SymEngine::hash_combine(seed, *p);
        return seed;
    }

    bool TwoElecEnergy::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const TwoElecEnergy>(o)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(o);
            return get_name() == op.get_name()
                && SymEngine::unified_eq(get_vec(), op.get_vec())
                && SymEngine::unified_eq(derivatives_, op.derivatives_)
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
                result = SymEngine::unified_compare(derivatives_, op.derivatives_);
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
    //    args.insert(args.end(), derivatives_.begin(), derivatives_.end());
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
        // tr(contr(g, inner_->diff(s)), outer_)
        auto diff_inner = SymEngine::rcp_dynamic_cast<const ElectronicState>(
            inner_->diff(s)
        );
        auto op_diff_inner = diff_inner->compare(*outer_)<=0
            ? SymEngine::make_rcp<const TwoElecEnergy>(
                  get_name(), diff_inner, outer_, dependencies_, derivatives_
              )
            : SymEngine::make_rcp<const TwoElecEnergy>(
                  get_name(), outer_, diff_inner, dependencies_, derivatives_
              );
        // tr(contr(g, inner_), outer_->diff(s))
        auto diff_outer = SymEngine::rcp_dynamic_cast<const ElectronicState>(
            outer_->diff(s)
        );
        auto op_diff_outer = inner_->compare(*diff_outer)<=0
            ? SymEngine::make_rcp<const TwoElecEnergy>(
                  get_name(), inner_, diff_outer, dependencies_, derivatives_
              )
            : SymEngine::make_rcp<const TwoElecEnergy>(
                  get_name(), diff_outer, inner_, dependencies_, derivatives_
              );
        auto op_diff_state = SymEngine::add(op_diff_inner, op_diff_outer);
        auto max_order = find_dependency(dependencies_, s);
        if (max_order > 0) {
            auto order = derivatives_.count(s) + 1;
            if (order <= max_order) {
                auto derivatives = derivatives_;
                derivatives.insert(s);
                return SymEngine::add(
                    // tr(contr(g->diff(s), inner_), outer_)
                    SymEngine::make_rcp<const TwoElecEnergy>(
                        get_name(),
                        inner_,
                        outer_,
                        dependencies_,
                        derivatives
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
