#include <symengine/symengine_assert.h>

#include "Tinned/OneElecDensity.hpp"

namespace Tinned
{
    OneElecDensity::OneElecDensity(
        const std::string& name,
        const SymEngine::multiset_basic& derivative
    ) : ElectronicState(name, derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::RCP<const SymEngine::Basic> OneElecDensity::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto derivative = derivative_;
        derivative.insert(s);
        return SymEngine::make_rcp<const OneElecDensity>(
            get_name(),
            derivative
        );
    }
}
