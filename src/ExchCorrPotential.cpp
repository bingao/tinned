#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/ExchCorrPotential.hpp"

namespace Tinned
{
    ExchCorrPotential::ExchCorrPotential(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::multiset_basic& derivative
    ) : SymEngine::MatrixSymbol(name),
        state_(state),
        derivative_(derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ExchCorrPotential::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine<const ElectronicState>(seed, *state_);
        for (auto& p: derivative_) {
            SymEngine::hash_combine<SymEngine::Basic>(seed, *p);
        }
        return seed;
    }

    bool ExchCorrPotential::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const ExchCorrPotential>(o)) {
                auto& op = SymEngine::down_cast<const ExchCorrPotential&>(o);
                // First check the electronic state
                if (not state_->__eq__(*op.state_)) return false;
                // Secondly check the derivatives
                return SymEngine::unified_eq(derivative_, op.derivative_);
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
                return SymEngine::unified_compare(derivative_, op.derivative_);
            }
            else {
                return result;
            }
        }
        return result;
    }

    SymEngine::vec_basic ExchCorrPotential::get_args() const
    {
        auto args = SymEngine::vec_basic({state_});
        args.insert(args.end(), derivative_.begin(), derivative_.end());
        return args;
    }

    SymEngine::RCP<const SymEngine::Basic> ExchCorrPotential::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto derivative = derivative_;
        derivative.insert(s);
        auto xc = SymEngine::make_rcp<const ExchCorrPotential>(
            SymEngine::MatrixSymbol::get_name(),
            state_,
            derivative
        );
        return xc;
    }
}
