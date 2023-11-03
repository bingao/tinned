#include <symengine/mul.h>
#include <symengine/symengine_casts.h>

#include "Tinned/ExchCorrEnergy.hpp"

namespace Tinned
{
    ExchCorrEnergy::ExchCorrEnergy(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::RCP<const OneElecOperator>& Omega,
        const SymEngine::RCP<const NonElecFunction>& weight,
        const unsigned int order
    ) : SymEngine::FunctionWrapper(name, SymEngine::vec_basic({weight, state, Omega})),
        energy_(SymEngine::mul(weight, make_exc_density(state, Omega, order)))
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    ExchCorrEnergy::ExchCorrEnergy(
        const ExchCorrEnergy& other,
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) : SymEngine::FunctionWrapper(other.get_name(), other.get_args()),
        energy_(other.energy_->diff(s))
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    ExchCorrEnergy::ExchCorrEnergy(
        const ExchCorrEnergy& other,
        const SymEngine::RCP<const SymEngine::Basic>& energy
    ) : SymEngine::FunctionWrapper(other.get_name(), other.get_args()),
        energy_(energy)
    {
        // XC energy or its derivatives must be either `SymEngine::Mul` or
        // `SymEngine::Add`
        SYMENGINE_ASSERT(
            SymEngine::is_a_sub<const SymEngine::Mul>(*energy) ||
            SymEngine::is_a_sub<const SymEngine::Add>(*energy)
        )
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ExchCorrEnergy::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::FunctionWrapper::__hash__();
        SymEngine::hash_combine(seed, *energy_);
        return seed;
    }

    bool ExchCorrEnergy::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::FunctionWrapper::__eq__(o)) {
            if (SymEngine::is_a_sub<const ExchCorrEnergy>(o)) {
                auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(o);
                return energy_->__eq__(*op.energy_);
            }
        }
        return false;
    }

    int ExchCorrEnergy::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ExchCorrEnergy>(o))
        int result = SymEngine::FunctionWrapper::compare(o);
        if (result == 0) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(o);
            return energy_->compare(*op.energy_);
        }
        return result;
    }

    SymEngine::RCP<const SymEngine::Basic> ExchCorrEnergy::create(
        const SymEngine::vec_basic &v
    ) const
    {
        throw SymEngine::NotImplementedError("ExchCorrEnergy::create() not implemented");
    }

    SymEngine::RCP<const SymEngine::Number> ExchCorrEnergy::eval(long bits) const
    {
        throw SymEngine::NotImplementedError("ExchCorrEnergy::eval() not implemented");
    }

    SymEngine::RCP<const SymEngine::Basic> ExchCorrEnergy::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        //return SymEngine::make_rcp<const ExchCorrEnergy>(
        //    *this,
        //    energy_->diff(s)
        //);
        return SymEngine::make_rcp<const ExchCorrEnergy>(*this, s);
    }
}
