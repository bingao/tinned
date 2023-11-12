#include <string>

#include <symengine/constants.h>

#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOverlap.hpp"

namespace Tinned
{
    TemporumOverlap::TemporumOverlap(const PertDependency& dependencies)
        : SymEngine::MatrixSymbol(std::string("T"))
    {
        // We use NonElecFunction to represent basis functions on bra and ket
        auto bra = SymEngine::make_rcp<const TemporumOperator>(
            SymEngine::make_rcp<const NonElecFunction>("<bra|", dependencies),
            TemporumType::Bra
        );
        auto ket = SymEngine::make_rcp<const TemporumOperator>(
            SymEngine::make_rcp<const NonElecFunction>("|ket>", dependencies),
            TemporumType::Ket
        );
        braket_ = SymEngine::matrix_mul({bra, ket});
        SYMENGINE_ASSIGN_TYPEID()
    }

    TemporumOverlap::TemporumOverlap(const SymEngine::RCP<const SymEngine::Basic>& braket)
        : SymEngine::MatrixSymbol(std::string("T")),
          braket_(braket)
    {
        SYMENGINE_ASSIGN_TYPEID()
        // If braket is an object of class SymEngine::MatrixAdd, each of its
        // terms should be an object of class SymEngine::MatrixMul, but we do
        // not bother to check them here
        SYMENGINE_ASSERT(
            SymEngine::is_a<const SymEngine::MatrixAdd>(*braket) ||
            SymEngine::is_a<const SymEngine::MatrixMul>(*braket)
        )
    }

    SymEngine::hash_t TemporumOverlap::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *braket_);
        return seed;
    }

    bool TemporumOverlap::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const TemporumOverlap>(o)) {
            auto& op = SymEngine::down_cast<const TemporumOverlap&>(o);
            return get_name() == op.get_name() && braket_->__eq__(*op.braket_);
        }
        return false;
    }

    int TemporumOverlap::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TemporumOverlap>(o))
        auto& op = SymEngine::down_cast<const TemporumOverlap&>(o);
        if (get_name() == op.get_name()) {
            return braket_->compare(*op.braket_);
        }
        else {
            return get_name() < op.get_name() ? -1 : 1;
        }
    }

    //SymEngine::vec_basic TemporumOverlap::get_args() const
    //{
    //    return SymEngine::vec_basic({braket_});
    //}

    SymEngine::RCP<const SymEngine::Basic> TemporumOverlap::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        auto result = braket_->diff(s);
        // `NonElecFunction` should return `SymEngine::zero`
        if (result->__eq__(*make_zero_operator()) || result->__eq__(*SymEngine::zero)) {
            return make_zero_operator();
        }
        else {
            return SymEngine::make_rcp<const TemporumOverlap>(result);
        }
    }
}
