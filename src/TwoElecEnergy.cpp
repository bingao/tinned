#include <symengine/add.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>
#include <symengine/matrices/matrix_add.h>

#include "Tinned/TwoElecEnergy.hpp"

namespace Tinned
{
    TwoElecEnergy::TwoElecEnergy(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& inner,
        const SymEngine::RCP<const ElectronicState>& outer,
        const PertDependency& dependencies,
        const SymEngine::multiset_basic& derivatives
    ) : SymEngine::FunctionWrapper(name, SymEngine::vec_basic({}))
    {
        // For numerical evaulation, one may put fewer density matrices in the
        // innermost loop. But we do not know the number of components of each
        // perturbation in Tinned, we also do not know the number of components
        // of high order differentiation with respect to a same perturbation.
        // So we simply rearrange inner_ and outer_ according to their
        // `compare` function, which helps the collection of equal
        // `TwoElecEnergy` objects.
        if (inner->compare(*outer)<=0) {
            G_ = SymEngine::make_rcp<const TwoElecOperator>(
                name, inner, dependencies, derivatives
            );
            outer_ = outer;
        }
        else {
            G_ = SymEngine::make_rcp<const TwoElecOperator>(
                name, outer, dependencies, derivatives
            );
            outer_ = inner;
        }
        SYMENGINE_ASSIGN_TYPEID()
    }

    TwoElecEnergy::TwoElecEnergy(
        const SymEngine::RCP<const TwoElecOperator>& G,
        const SymEngine::RCP<const ElectronicState>& outer
    ) : SymEngine::FunctionWrapper(G->get_name(), SymEngine::vec_basic({}))
    {
        auto inner = G->get_state();
        if (inner->compare(*outer)<=0) {
            G_ = G;
            outer_ = outer;
        }
        else {
            G_ = SymEngine::make_rcp<const TwoElecOperator>(
                G->get_name(), outer, G->get_dependencies(), G->get_derivatives()
            );
            outer_ = inner;
        }
    }

    SymEngine::hash_t TwoElecEnergy::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::FunctionWrapper::__hash__();
        SymEngine::hash_combine(seed, *G_);
        SymEngine::hash_combine(seed, *outer_);
        return seed;
    }

    bool TwoElecEnergy::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const TwoElecEnergy>(o)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(o);
            return get_name()==op.get_name()
                && G_->__eq__(*op.G_) && outer_->__eq__(*op.outer_);
        }
        return false;
    }

    int TwoElecEnergy::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TwoElecEnergy>(o))
        auto& op = SymEngine::down_cast<const TwoElecEnergy&>(o);
        if (get_name()==op.get_name()) {
            int result = G_->compare(*op.G_);
            return result==0 ? outer_->compare(*op.outer_) : result;
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic TwoElecEnergy::get_args() const
    {
        return SymEngine::vec_basic({G_, outer_});
    }

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
        // tr(contr(g_, inner_), outer_->diff(s))
        auto terms = SymEngine::vec_basic({
            SymEngine::make_rcp<const TwoElecEnergy>(
                G_, SymEngine::rcp_dynamic_cast<const ElectronicState>(outer_->diff(s))
            )
        });
        // tr(contr(g_, inner_)->diff(s), outer_), where contr(g_, inner_)->diff(s)
        // is either `TwoElecOperator` or `SymEngine::MatrixAdd`, see
        // `src/TwoElecOperator.cpp`
        auto diff_G = G_->diff(s);
        if (SymEngine::is_a_sub<const TwoElecOperator>(*diff_G)) {
            terms.push_back(SymEngine::make_rcp<const TwoElecEnergy>(
                SymEngine::rcp_dynamic_cast<const TwoElecOperator>(diff_G), outer_
            ));
        }
        else {
            SYMENGINE_ASSERT(SymEngine::is_a_sub<const SymEngine::MatrixAdd>(*diff_G))
            auto op = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixAdd>(diff_G);
            for (const auto& arg: op->get_args()) {
                SYMENGINE_ASSERT(SymEngine::is_a_sub<const TwoElecOperator>(*arg))
                terms.push_back(SymEngine::make_rcp<const TwoElecEnergy>(
                    SymEngine::rcp_dynamic_cast<const TwoElecOperator>(arg), outer_
                ));
            }
        }
        return SymEngine::add(terms);
    }
}
