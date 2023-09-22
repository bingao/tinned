#include <string>

#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/Perturbation.hpp"

namespace Tinned
{
    Perturbation::Perturbation(
        const char name,
        const std::size_t dimension
    ) : SymEngine::Symbol(std::string(1, name)), dimension_(dimension)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    // We use ASCII code as the hash of a perturbation by noticing that the
    // name of the perturbation is a single character
    SymEngine::hash_t Perturbation::__hash__() const
    {
        return SymEngine::hash_t(get_name()[0]);
    }

    bool Perturbation::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::Symbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const Perturbation>(o))
                return dimension_ == SymEngine::down_cast<const Perturbation &>(o).dimension_;
        }
        return false;
    }

    int Perturbation::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(o))
        int result = SymEngine::Symbol::compare(o);
        if (result == 0) {
            //const Perturbation &s = SymEngine::down_cast<const Perturbation &>(o);
            const Perturbation& op = SymEngine::down_cast<const Perturbation &>(o);
            return dimension_ < op.dimension_ ? -1 : 1;
        }
        return result;
    }
}
