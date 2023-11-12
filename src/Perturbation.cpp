#include <iostream>

#include <symengine/complex.h>
//#include <symengine/complex_double.h>
//#include <symengine/complex_mpc.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/Perturbation.hpp"

namespace Tinned
{
    Perturbation::Perturbation(
        const std::string& name,
        const SymEngine::RCP<const SymEngine::Number>& frequency,
        const std::set<std::size_t>& components
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
        for (auto& c: components_) SymEngine::hash_combine(seed, c);
        return seed;
    }

    bool Perturbation::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const Perturbation>(o)) {
            auto& s = SymEngine::down_cast<const Perturbation&>(o);
            return get_name() == s.get_name()
                && frequency_->__eq__(*s.frequency_)
                && SymEngine::ordered_eq(components_, s.components_);
        }
        return false;
    }

    int Perturbation::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const Perturbation>(o))
        auto& s = SymEngine::down_cast<const Perturbation&>(o);
        if (get_name() == s.get_name()) {
            int result;
            // Some subclasses of SymEngine::Number cannot be compared
            // directly, so we take their difference and compare
            auto diff = SymEngine::subnum(frequency_, s.frequency_);
            if (diff->is_complex()) {
                auto diff_cmplx = SymEngine::rcp_dynamic_cast<const SymEngine::ComplexBase>(diff);
                auto diff_real = diff_cmplx->real_part();
                if (diff_real->is_zero()) {
                    auto diff_imag = diff_cmplx->imaginary_part();
                    if (diff_imag->is_zero()) {
                        result = 0;
                    }
                    else if (diff_imag->is_negative()) {
                        result = -1;
                    }
                    else {
                        result = 1;
                    }
                }
                else if (diff_real->is_negative()) {
                    result = -1;
                }
                else {
                    result = 1;
                }
            }
            else {
                if (diff->is_zero()) {
                    result = 0;
                }
                else if (diff->is_negative()) {
                    result = -1;
                }
                else {
                    result = 1;
                }
            }
            return result == 0
                ? SymEngine::ordered_compare(components_, s.components_)
                : result;
        }
        else {
            return get_name() < s.get_name() ? -1 : 1;
        }
    }
}
