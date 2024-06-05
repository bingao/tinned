#include <string>

//#include <symengine/integer.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>

#include "Tinned/TemporumOperator.hpp"
#include "Tinned/ZeroOperator.hpp"

namespace Tinned
{
    TemporumOperator::TemporumOperator(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& target,
        const TemporumType type
    ) : SymEngine::MatrixSymbol(
            type==TemporumType::Ket ? std::string("i*dt") : std::string("-i*dt")
        ),
        target_(target),
        type_(type)
    {
        SYMENGINE_ASSIGN_TYPEID()
        // Currently we allow only the following classes as the target that the
        // time differentiation operator acts on
        SYMENGINE_ASSERT(
            SymEngine::is_a_sub<const ElectronicState>(*target) ||
            SymEngine::is_a_sub<const OneElecOperator>(*target)
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
        if (SymEngine::is_a_sub<const TemporumOperator>(o)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(o);
            return get_name()==op.get_name()
                && type_==op.type_
                && target_->__eq__(*op.target_);
        }
        return false;
    }

    int TemporumOperator::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TemporumOperator>(o))
        auto& op = SymEngine::down_cast<const TemporumOperator&>(o);
        if (get_name()==op.get_name()) {
            int result = SymEngine::unified_compare(
                static_cast<int>(type_),
                static_cast<int>(op.type_)
            );
            return result==0 ? target_->compare(*op.target_) : result;
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic TemporumOperator::get_args() const
    {
        return SymEngine::vec_basic({
            //SymEngine::integer(static_cast<int>(type_)),
            target_
        });
    }

    SymEngine::RCP<const SymEngine::Basic> TemporumOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto result = target_->diff(s);
        if (result->__eq__(*make_zero_operator())) {
            return make_zero_operator();
        }
        else {
            return SymEngine::make_rcp<const TemporumOperator>(
                SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(result), type_
            );
        }
    }
}
