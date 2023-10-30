#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/matrices/matrix_mul.h>

#include "Tinned/ExchCorrPotential.hpp"

namespace Tinned
{
    ExchCorrPotential::ExchCorrPotential(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const NonElecFunction>& weight
    ) : SymEngine::MatrixSymbol(name),
        state_(state),
        Omega_(Omega),
        weight_(weight),
        potential_(SymEngine::matrix_mul(
            SymEngine::vec_basic({make_vxc(name, state, Omega, weight), Omega})
        ))
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    ExchCorrPotential::ExchCorrPotential(
        const ExchCorrPotential& other,
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) : SymEngine::MatrixSymbol(other.get_name()),
        state_(other.state_),
        Omega_(other.Omega_),
        weight_(other.weight_),
        potential_(other.potential_->diff(s))
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    ExchCorrPotential::ExchCorrPotential(
        const ExchCorrPotential& other,
        const SymEngine::RCP<const SymEngine::Basic>& potential
    ) : SymEngine::MatrixSymbol(other.get_name()),
        state_(other.state_),
        Omega_(other.Omega_),
        weight_(other.weight_),
        potential_(potential)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ExchCorrPotential::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *state_);
        SymEngine::hash_combine(seed, *weight_);
        SymEngine::hash_combine(seed, *Omega_);
        SymEngine::hash_combine(seed, *potential_);
        return seed;
    }

    bool ExchCorrPotential::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const ExchCorrPotential>(o)) {
                auto& op = SymEngine::down_cast<const ExchCorrPotential&>(o);
                return state_->__eq__(*op.state_)
                    && weight_->__eq__(*op.weight_)
                    && Omega_->__eq__(*op.Omega_)
                    && potential_->__eq__(*op.potential_);
            }
        }
        return false;
    }

    int ExchCorrPotential::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ExchCorrPotential>(o))
        int result = SymEngine::MatrixSymbol::compare(o);
        if (result == 0) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(o);
            result = state_->compare(*op.state_);
            if (result == 0) {
                result = weight_->compare(*op.weight_);
                if (result == 0) {
                    result = Omega_->compare(*op.Omega_);
                    return result == 0 ? potential_->compare(*op.potential_) : result;
                }
                else {
                    return result;
                }
            }
            else {
                return result;
            }
        }
        return result;
    }

    SymEngine::vec_basic ExchCorrPotential::get_args() const
    {
        return SymEngine::vec_basic({weight_, state_, Omega_});
    }

    SymEngine::RCP<const SymEngine::Basic> ExchCorrPotential::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        return SymEngine::make_rcp<const ExchCorrPotential>(*this, s);
    }
}
