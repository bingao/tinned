#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/ElectronicState.hpp"

namespace Tinned
{
    ElectronicState::ElectronicState(
        const std::string& name,
        const SymEngine::multiset_basic& derivatives
    ) : SymEngine::MatrixSymbol(name),
        derivatives_(derivatives)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ElectronicState::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        for (auto& p: derivatives_) SymEngine::hash_combine(seed, *p);
        return seed;
    }

    bool ElectronicState::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const ElectronicState>(o)) {
            auto& state = SymEngine::down_cast<const ElectronicState&>(o);
            return get_name() == state.get_name()
                && SymEngine::unified_eq(derivatives_, state.derivatives_);
        }
        return false;
    }

    int ElectronicState::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ElectronicState>(o))
        auto& state = SymEngine::down_cast<const ElectronicState&>(o);
        if (get_name() == state.get_name()) {
            return SymEngine::unified_compare(derivatives_, state.derivatives_);
        }
        else {
            return get_name() < state.get_name() ? -1 : 1;
        }
    }

    //SymEngine::vec_basic ElectronicState::get_args() const
    //{
    //    if (derivatives_.empty()) {
    //        return {};
    //    }
    //    else {
    //        return SymEngine::vec_basic(derivatives_.begin(), derivatives_.end());
    //    }
    //}

    //SymEngine::RCP<const SymEngine::Basic> ElectronicState::diff_impl(
    //    const SymEngine::RCP<const SymEngine::Symbol>& s
    //) const
    //{
    //    auto derivatives = derivatives_;
    //    derivatives.insert(s);
    //    return SymEngine::make_rcp<const ElectronicState>(
    //        get_name(),
    //        derivatives
    //    );
    //}
}
