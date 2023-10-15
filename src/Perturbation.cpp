#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/Perturbation.hpp"

namespace Tinned
{
    Perturbation::Perturbation(
        const std::string& name,
        const SymEngine::RCP<const SymEngine::Number>& frequency,
        const std::set<std::size_t> components
    ) : SymEngine::Symbol(name),
        frequency_(frequency),
        components_(components)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t Perturbation::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::Symbol::__hash__();
        SymEngine::hash_combine(seed, *frequency_);
        for (auto& c: components_) SymEngine::hash_combine<std::size_t>(seed, c);
        return seed;
    }

    bool Perturbation::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::Symbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const Perturbation>(o)) {
                auto& s = SymEngine::down_cast<const Perturbation&>(o);
                return frequency_->__eq__(*s.frequency_)
                    && SymEngine::ordered_eq(components_, s.components_);
            }
        }
        return false;
    }

    int Perturbation::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(o))
        int result = SymEngine::Symbol::compare(o);
        if (result == 0) {
            auto& s = SymEngine::down_cast<const Perturbation&>(o);
            result = frequency_->compare(*s.frequency_);
            return result == 0
                ? SymEngine::ordered_compare(components_, s.components_)
                : result;
        }
        return result;
    }
}
