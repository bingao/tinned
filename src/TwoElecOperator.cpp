#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/matrices/matrix_add.h>

#include "Tinned/TwoElecOperator.hpp"

namespace Tinned
{
    TwoElecOperator::TwoElecOperator(
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

    SymEngine::hash_t TwoElecOperator::__hash__() const
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

    bool TwoElecOperator::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::MatrixSymbol::__eq__(o)) {
            if (SymEngine::is_a_sub<const TwoElecOperator>(o)) {
                const TwoElecOperator& op = SymEngine::down_cast<const TwoElecOperator &>(o);
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

    int TwoElecOperator::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const TwoElecOperator>(o))
        int result = SymEngine::MatrixSymbol::compare(o);
        if (result == 0) {
            const TwoElecOperator& op = SymEngine::down_cast<const TwoElecOperator &>(o);
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

    SymEngine::vec_basic TwoElecOperator::get_args() const
    {
        auto args = SymEngine::vec_basic({state_});
        auto deps = to_vec_basic(dependencies_);
        args.insert(args.end(), deps.begin(), deps.end());
        args.insert(args.end(), derivative_.begin(), derivative_.end());
        return args;
    }

    SymEngine::RCP<const SymEngine::MatrixExpr> TwoElecOperator::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        // contr(g, D->diff(s))
        auto diff_state = SymEngine::rcp_dynamic_cast<const ElectronState>(
            state_->diff(s)
        );
        auto op_diff_state = SymEngine::make_rcp<const TwoElecOperator>(
            SymEngine::MatrixSymbol::get_name(),
            diff_state,
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
                    SymEngine::make_rcp<const TwoElecOperator>(
                        SymEngine::MatrixSymbol::get_name(),
                        state_,
                        dependencies_,
                        derivative
                    ),
                    op_diff_state
                });
                return op;
            }
            else {
                return op_diff_state;
            }
        }
        else {
            return op_diff_state;
        }
    }
}
