#include <string>

#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/integer.h>

#include "Tinned/TemporumOperator.hpp"

namespace Tinned
{
    TemporumOperator::TemporumOperator(
        const SymEngine::RCP<const Basic>& target,
        const TemporumType type
    ) : SymEngine::MatrixSymbol(
            type == TemporumType::Ket ? std::string("i*dt") : std::string("-i*dt")
        ),
        target_(target),
        type_(type)
    {
        SYMENGINE_ASSIGN_TYPEID()
        // Currently we only allow the following classes as the target that the
        // time differentiation operator acts on
        SYMENGINE_ASSERT(
            SymEngine::is_a_sub<const ElectronicState>(*target) ||
            SymEngine::is_a_sub<const OneElecOperator>(*target) ||
            SymEngine::is_a_sub<const NonElecFunction>(*target)
        )
    }

    SymEngine::hash_t TemporumOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *target_);
        SymEngine::hash_combine(seed, type_);
        return seed;
    }

    bool TemporumOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const TemporumOperator>(o)) {
                auto& op = SymEngine::down_cast<const TemporumOperator&>(o);
                return type_ == op.type_ && target_->__eq__(*op.target_);
            }
        }
        return false;
    }

    int TemporumOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TemporumOperator>(o))
        int result = SymEngine::MatrixSymbol::compare(o);
        if (result == 0) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(o);
            result = SymEngine::unified_compare(
                static_cast<int>(type_),
                static_cast<int>(op.type_)
            );
            return result == 0 ? target_->compare(*op.target_) : result;
        }
        return result;
    }

    SymEngine::vec_basic TemporumOperator::get_args() const
    {
        return SymEngine::vec_basic({
            SymEngine::integer(static_cast<int>(type_)),
            target_
        });
    }

    SymEngine::RCP<const SymEngine::Basic> TemporumOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto result = target_->diff(s);
        if (zero_operator()->__eq__(*result)) {
            return zero_operator();
        }
        else {
            return SymEngine::make_rcp<const TemporumOperator>(result, type_);
        }
    }
}
