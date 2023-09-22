#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/matrices/matrix_add.h>

#include "Tinned/TwoElectronOperator.hpp"

namespace Tinned
{
    TwoElectronOperator::TwoElectronOperator(
        const std::string& name,
        const SymEngine::RCP<const ElectronState>& state,
        const PertDependency& dependencies,
        const SymEngine::multiset_basic& derivative
    ) : SymEngine::MatrixSymbol(name),
        state_(state),
        dependencies_(dependencies),
        derivative_(derivative)
    {
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t TwoElectronOperator::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine<const ElectronState>(seed, *state_);
        for (auto& dep: dependencies_) {
            SymEngine::hash_combine<const Perturbation>(seed, *dep.first);
            SymEngine::hash_combine<unsigned int>(seed, dep.second);
        }
        for (auto& p: derivative_) {
            SymEngine::hash_combine<SymEngine::Basic>(seed, *p);
        }
        return seed;
    }

    bool TwoElectronOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const TwoElectronOperator>(o)) {
                const TwoElectronOperator& op = SymEngine::down_cast<const TwoElectronOperator &>(o);
                // First check the electron state
                if (not state_->__eq__(*op.state_)) return false;
                // Secondly check the derivatives
                if (not SymEngine::unified_eq(derivative_, op.derivative_)) return false;
                // Thirdly we check the perturbation dependencies
                return eq_dependency(dependencies_, op.dependencies_);
            }
        }
        return false;
    }

    int TwoElectronOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TwoElectronOperator>(o))
        int result = SymEngine::MatrixSymbol::compare(o);
        if (result == 0) {
            const TwoElectronOperator& op = SymEngine::down_cast<const TwoElectronOperator &>(o);
            result = state_->compare(*op.state_);
            if (result == 0) {
                result = SymEngine::unified_compare(derivative_, op.derivative_);
                if (result == 0) {
                    return SymEngine::ordered_compare(dependencies_, op.dependencies_);
                }
                else {
                    return result;
                }
            }
            else {
                return result;
            }
        }
        return result;
    }

    SymEngine::vec_basic TwoElectronOperator::get_args() const
    {
        auto args = SymEngine::vec_basic({state_});
        auto deps = to_vec_basic(dependencies_);
        args.insert(args.end(), deps.begin(), deps.end());
        args.insert(args.end(), derivative_.begin(), derivative_.end());
        return args;
    }

    SymEngine::RCP<const SymEngine::MatrixExpr> TwoElectronOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        // contr(g, D->diff(s))
        auto op_ds = SymEngine::make_rcp<const TwoElectronOperator>(
            SymEngine::MatrixSymbol::get_name(),
            state_->diff(s),
            dependencies_,
            derivative_
        );
        auto max_order = find_dependency(dependencies_, s);
        if (max_order > 0) {
            auto order = derivative_.count(s) + 1;
            if (order <= max_order) {
                // Return contr(g->diff(s), D) + contr(g, D->diff(s))
                auto derivative = derivative_;
                derivative.insert(s);
                auto op = SymEngine::matrix_add({
                    SymEngine::make_rcp<const TwoElectronOperator>(
                        SymEngine::MatrixSymbol::get_name(),
                        state_,
                        dependencies_,
                        derivative
                    ),
                    op_ds
                });
                return op;
            }
            else {
                return op_ds;
            }
        }
        else {
            return op_ds;
        }
    }
}
