#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/ExchCorrPotential.hpp"

namespace Tinned
{
    ExchCorrPotential::ExchCorrPotential(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const PertDependency& dependencies,
        const SymEngine::RCP<const NonElecFunction>& weight
    ) : SymEngine::MatrixSymbol(name),
        potential_(SymEngine::matrix_mul(
            SymEngine::vec_basic({
                // Exchange-correlation potential, or the functional derivative
                // of exchange-correlation energy
                SymEngine::make_rcp<const ExchCorrEnergy>(
                    std::string("vxc"),
                    state,
                    weight,
                    1
                ),
                // Overlap distribution
                SymEngine::make_rcp<const OneElecOperator>(
                    std::string("Omega"),
                    dependencies
                )
            })
        ))
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    ExchCorrPotential::ExchCorrPotential(
        const ExchCorrPotential& other,
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) : SymEngine::MatrixSymbol(other.get_name()),
        potential_(other.potential_->diff(s))
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ExchCorrPotential::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *potential_);
        return seed;
    }

    bool ExchCorrPotential::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const ExchCorrPotential>(o)) {
                auto& op = SymEngine::down_cast<const ExchCorrPotential&>(o);
                return potential_->__eq__(*op.potential_);
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
            return potential_->compare(*op.potential_);
        }
        return result;
    }

    SymEngine::vec_basic ExchCorrPotential::get_args() const
    {
        return SymEngine::vec_basic({state_, weight_, Omega_});
    }

    SymEngine::RCP<const SymEngine::Basic> ExchCorrPotential::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        return SymEngine::make_rcp<const ExchCorrPotential>(*this, s);
    }
}
