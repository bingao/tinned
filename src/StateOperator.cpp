#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/StateOperator.hpp"

namespace Tinned
{
    StateOperator::StateOperator(
        const std::string& name,
        const SymEngine::RCP<const StateVector>& state
    ) : SymEngine::MatrixSymbol(name),
        state_(state)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t StateOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *state_);
        return seed;
    }

    bool StateOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const StateOperator>(o)) {
            auto& op = SymEngine::down_cast<const StateOperator&>(o);
            return get_name() == op.get_name() && state_->__eq__(*op.state_);
        }
        return false;
    }

    int StateOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const StateOperator>(o))
        auto& op = SymEngine::down_cast<const StateOperator&>(o);
        if (get_name() == op.get_name()) {
            return state_->compare(*op.state_);
        }
        else {
            return get_name() < op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic StateOperator::get_args() const
    {
        return SymEngine::vec_basic({state_});
    }

    SymEngine::RCP<const SymEngine::Basic> StateOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto diff_state = SymEngine::rcp_dynamic_cast<const StateVector>(
            state_->diff(s)
        );
        if (diff_state->__eq__(*make_zero_operator())) {
            return make_zero_operator();
        }
        else {
            return SymEngine::make_rcp<const StateOperator>(get_name(), diff_state);
        }
    }
}
