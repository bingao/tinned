#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/ElectronicState.hpp"

namespace Tinned
{
    ElectronicState::ElectronicState(
        const std::string& name,
        const SymEngine::multiset_basic& derivative
    ) : SymEngine::MatrixSymbol(name),
        derivative_(derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ElectronicState::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        for (auto& p: derivative_) {
            SymEngine::hash_combine<SymEngine::Basic>(seed, *p);
        }
        return seed;
    }

    bool ElectronicState::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const ElectronicState>(o)) {
                auto& state = SymEngine::down_cast<const ElectronicState&>(o);
                // Check the derivatives
                return SymEngine::unified_eq(derivative_, state.derivative_);
            }
        }
        return false;
    }

    int ElectronicState::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ElectronicState>(o))
        int result = SymEngine::MatrixSymbol::compare(o);
        if (result == 0) {
            auto& state = SymEngine::down_cast<const ElectronicState&>(o);
            return SymEngine::unified_compare(derivative_, state.derivative_);
        }
        return result;
    }

    SymEngine::vec_basic ElectronicState::get_args() const
    {
        if (derivative_.empty()) {
            return {};
        }
        else {
            return SymEngine::vec_basic(derivative_.begin(), derivative_.end());
        }
    }

    //SymEngine::RCP<const SymEngine::Basic> ElectronicState::diff_impl(
    //    const SymEngine::RCP<const SymEngine::Symbol>& s
    //) const
    //{
    //    auto derivative = derivative_;
    //    derivative.insert(s);
    //    auto state = SymEngine::make_rcp<const ElectronicState>(
    //        SymEngine::MatrixSymbol::get_name(),
    //        derivative
    //    );
    //    return state;
    //}
}
