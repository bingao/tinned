#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/PerturbedParameter.hpp"

namespace Tinned
{
    PerturbedParameter::PerturbedParameter(
        const std::string& name,
        const SymEngine::multiset_basic& derivatives
    ) : SymEngine::MatrixSymbol(name),
        derivatives_(derivatives)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t PerturbedParameter::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        for (auto& p: derivatives_) SymEngine::hash_combine(seed, *p);
        return seed;
    }

    bool PerturbedParameter::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const PerturbedParameter>(o)) {
            auto& op = SymEngine::down_cast<const PerturbedParameter&>(o);
            return get_name()==op.get_name()
                && SymEngine::unified_eq(derivatives_, op.derivatives_);
        }
        return false;
    }

    int PerturbedParameter::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const PerturbedParameter>(o))
        auto& op = SymEngine::down_cast<const PerturbedParameter&>(o);
        if (get_name()==op.get_name()) {
            return SymEngine::unified_compare(derivatives_, op.derivatives_);
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::RCP<const SymEngine::Basic> PerturbedParameter::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto derivatives = derivatives_;
        derivatives.insert(s);
        return SymEngine::make_rcp<const PerturbedParameter>(
            get_name(),
            derivatives
        );
    }
}
