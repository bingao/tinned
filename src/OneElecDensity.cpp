#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/OneElecDensity.hpp"

namespace Tinned
{
    OneElecDensity::OneElecDensity(
        const std::string& name,
        const SymEngine::multiset_basic& derivative
    ) : ElectronState(name, derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::RCP<const SymEngine::MatrixExpr> OneElecDensity::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto derivative = derivative_;
        derivative.insert(s);
        auto state = SymEngine::make_rcp<const OneElecDensity>(
            SymEngine::MatrixSymbol::get_name(),
            derivative
        );
        return state;
    }
}
