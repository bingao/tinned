#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/ExchCorrEnergy.hpp"

namespace Tinned
{
    ExchCorrEnergy::ExchCorrEnergy(
        const std::string& name,
        const SymEngine::RCP<const ElectronicState>& state,
        const SymEngine::multiset_basic& derivative
    ) : SymEngine::FunctionWrapper(name, state),
        derivative_(derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ExchCorrEnergy::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::FunctionWrapper::__hash__();
        for (auto& p: derivative_) {
            SymEngine::hash_combine<SymEngine::Basic>(seed, *p);
        }
        return seed;
    }

    bool ExchCorrEnergy::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::FunctionWrapper::__eq__(o)) {
            if (SymEngine::is_a_sub<const ExchCorrEnergy>(o)) {
                auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(o);
                return SymEngine::unified_eq(derivative_, op.derivative_);
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
            return SymEngine::unified_compare(derivative_, op.derivative_);
        }
        return result;
    }

    SymEngine::vec_basic ExchCorrEnergy::get_args() const
    {
        auto args = SymEngine::FunctionWrapper::get_args();
        args.insert(args.end(), derivative_.begin(), derivative_.end());
        return args;
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
        auto derivative = derivative_;
        derivative.insert(s);
        auto xc = SymEngine::make_rcp<const ExchCorrEnergy>(
            SymEngine::FunctionWrapper::get_name(),
            get_state(),
            derivative
        );
        return xc;
    }
}
