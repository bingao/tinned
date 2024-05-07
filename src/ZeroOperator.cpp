#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/ZeroOperator.hpp"

namespace Tinned
{
    SymEngine::hash_t ZeroOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        return seed;
    }

    bool ZeroOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const ZeroOperator>(o)) {
            auto& op = SymEngine::down_cast<const ZeroOperator&>(o);
            return get_name()==op.get_name();
        }
        return false;
    }

    int ZeroOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ZeroOperator>(o))
        auto& op = SymEngine::down_cast<const ZeroOperator&>(o);
        if (get_name()==op.get_name()) {
            return 0;
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::RCP<const SymEngine::Basic> ZeroOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        return make_zero_operator();
    }
}
