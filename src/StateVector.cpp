#include <symengine/symengine_assert.h>

#include "Tinned/StateVector.hpp"

namespace Tinned
{
    StateVector::StateVector(
        const std::string& name,
        const SymEngine::multiset_basic& derivatives
    ) : ElectronicState(name, derivatives)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::RCP<const SymEngine::Basic> StateVector::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto derivatives = derivatives_;
        derivatives.insert(s);
        return SymEngine::make_rcp<const StateVector>(
            get_name(),
            derivatives
        );
    }
}
