#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/ZeroOperator.hpp"

namespace Tinned
{
    SymEngine::hash_t ZeroOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::ZeroMatrix::__hash__();
        return seed;
    }

    bool ZeroOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const ZeroOperator>(o)) {
            return true;
            //auto& op = SymEngine::down_cast<const ZeroOperator&>(o);
            //return get_name()==op.get_name();
        }
        return false;
    }

    int ZeroOperator::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ZeroOperator>(o))
        return 0;
        //auto& op = SymEngine::down_cast<const ZeroOperator&>(o);
        //if (get_name()==op.get_name()) {
        //    return 0;
        //}
        //else {
        //    return get_name()<op.get_name() ? -1 : 1;
        //}
    }
}
