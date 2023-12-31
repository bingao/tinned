#include <symengine/symengine_assert.h>

#include "Tinned/OneElecDensity.hpp"

namespace Tinned
{
    OneElecDensity::OneElecDensity(
        const std::string& name,
        const SymEngine::multiset_basic& derivatives
    ) : ElectronicState(name, derivatives)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::RCP<const SymEngine::Basic> OneElecDensity::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto derivatives = derivatives_;
        derivatives.insert(s);
        return SymEngine::make_rcp<const OneElecDensity>(
            get_name(),
            derivatives
        );
    }
}
