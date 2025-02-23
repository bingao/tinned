#include <symengine/integer.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/operators/CompositeFunction.hpp"

namespace Tinned
{
    CompositeFunction::CompositeFunction(
        const std::string& name,
        const SymEngine::RCP<const SymEngine::Basic> inner,
        const unsigned int order
    ) : SymEngine::FunctionWrapper(name, SymEngine::vec_basic({})),
        order_(order),
        inner_(inner)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t CompositeFunction::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::FunctionWrapper::__hash__();
        SymEngine::hash_combine(seed, order_);
        SymEngine::hash_combine(seed, *inner_);
        return seed;
    }

    bool CompositeFunction::__eq__(const SymEngine::Basic& o) const
    {
        // We did not introduce `type_code_id`'s for classes derived from
        // `FunctionWrapper`, such as `CompositeFunction` and
        // `NonElecFunction`, so `SymEngine::is_a<T>` will always be true for
        // the comparison of objects from these derived classes. The use of
        // `SymEngine::is_a_sub<T>` solves this problem, but it can not
        // distinguish, for example, class `CompositeFunction` and its derived
        // ones. So derived classes from those implemented in Tinned are not
        // recommended.
        if (SymEngine::is_a_sub<const CompositeFunction>(o)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(o);
            return get_name()==op.get_name()
                && order_==op.order_
                && inner_->__eq__(*op.inner_);
        }
        return false;
    }

    int CompositeFunction::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const CompositeFunction>(o))
        auto& op = SymEngine::down_cast<const CompositeFunction&>(o);
        if (get_name()==op.get_name()) {
            int result = SymEngine::unified_compare(order_, op.order_);
            return result==0 ? inner_->compare(*op.inner_) : result;
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic CompositeFunction::get_args() const
    {
        return SymEngine::vec_basic({SymEngine::integer(order_), inner_});
    }

    SymEngine::RCP<const SymEngine::Basic> CompositeFunction::create(
        const SymEngine::vec_basic &v
    ) const
    {
        throw SymEngine::NotImplementedError("CompositeFunction::create() not implemented");
    }

    SymEngine::RCP<const SymEngine::Number> CompositeFunction::eval(long bits) const
    {
        throw SymEngine::NotImplementedError("CompositeFunction::eval() not implemented");
    }

    SymEngine::RCP<const SymEngine::Basic> CompositeFunction::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        // Zero will be returned if the derivatives of the inner function
        // become zero
        return SymEngine::mul(
            SymEngine::make_rcp<const CompositeFunction>(
                get_name(), inner_, order_+1
            ),
            inner_->diff(s)
        );
    }
}
