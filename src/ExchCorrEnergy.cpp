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
    ) : SymEngine::FunctionWrapper(other.get_name(), other.get_args())
    {
        auto state = other.get_state();
        auto Omega = other.get_overlap_distribution();
        SymEngine::vec_basic terms = {};
        auto energy_map = extract_energy_map(energy);
        // `weight_map.first` is the (un)perturbed grid weight, and
        // `weight_map.second` is type `ExcDensityContractionMap`
        for (const auto& weight_map: energy_map) {
            // `exc_map.first` is the order of XC energy functional derivative
            // vector(s), and `exc_map.second` is the generalized density
            // vector(s)
            for (const auto& exc_map: weight_map.second) {
                if (exc_map.second.is_null()) {
                    terms.push_back(SymEngine::mul(
                        weight_map.first,
                        make_exc_density(state, Omega, exc_map.first)
                    ));
                }
                else {
                    terms.push_back(SymEngine::mul(SymEngine::vec_basic({
                        weight_map.first,
                        make_exc_density(state, Omega, exc_map.first),
                        exc_map.second
                    })));
                }
            }
        }
        SYMENGINE_ASSERT(!terms.empty())
        if (terms.size() == 1) {
            energy_ = terms[0];
        }
        else {
            energy_ = SymEngine::add(terms);
        }
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
        // Be careful of using base class' `__eq__()` method.
        // `SymEngine::FunctionSymbol::__eq__()` requires
        // `SymEngine::is_a<FunctionSymbol>(o)` which is not true here.
        if (SymEngine::is_a_sub<const ExchCorrEnergy>(o)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(o);
            return get_name() == op.get_name()
                && SymEngine::unified_eq(get_vec(), op.get_vec())
                && energy_->__eq__(*op.energy_);
        }
        return false;
    }

    int ExchCorrEnergy::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ExchCorrEnergy>(o))
        auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(o);
        if (get_name() == op.get_name()) {
            int result = SymEngine::unified_compare(get_vec(), op.get_vec());
            return result == 0 ? energy_->compare(*op.energy_) : result;
        }
        else {
            return get_name() < op.get_name() ? -1 : 1;
        }
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
