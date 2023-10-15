#include <symengine/symengine_exception.h>

#include "Tinned/RemoveVisitor.hpp"

namespace Tinned
{
    RemoveVisitor::RemoveVisitor(const SymEngine::vec_basic& symbols, bool equivalence)
        : symbols_(symbols)
    {
        if (equivalence) {
            to_remove_ = &eq_check;
        }
        else {
            to_remove_ = &neq_check;
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::Basic& x)
    {
        throw SymEngine::NotImplementedError(
            "RemoveVisitor::bvisit() not implemented for " + x.__str__()
        );
    }

    // Each upper level check if a zero matrix/symbol returned from lower level

    void RemoveVisitor::bvisit(const SymEngine::Symbol& x)
    {
        if (SymEngine::is_a_sub<const Perturbation>(x)) {
            if (to_remove_(x)) {
                result_ = SymEngine::RCP();
            }
            else {
                result_ = x.rcp_from_this();
            }
        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveVisitor::bvisit() not implemented for " + x.__str__()
            );
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        // We don't allow for the removal of derivative symbols, but only check
        // if the `NonElecFunction` can be removed as a whole
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            if (to_remove_(x)) {
                result_ = SymEngine::RCP();
            }
            else {
                result_ = x.rcp_from_this();
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            str_ = to_string(op.get_name(), op.get_state(), op.get_derivative());
        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveVisitor::bvisit() not implemented for " + x.__str__()
            );
        }
    }

    void RemoveVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            auto derivative = SymEngine::down_cast<const OneElecDensity&>(x).get_derivative();
            if (derivative.empty()) {
                str_ = SymEngine::down_cast<const OneElecDensity&>(x).get_name();
            }
            else {
                str_ = to_string(
                    SymEngine::down_cast<const OneElecDensity&>(x).get_name(),
                    derivative
                );
            }
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            auto name = SymEngine::down_cast<const OneElecOperator&>(x).get_name();
            auto derivative = SymEngine::down_cast<const OneElecOperator&>(x).get_derivative();
            if (derivative.empty()) {
                str_ = to_string(
                    name,
                    SymEngine::down_cast<const OneElecOperator&>(x).get_dependencies()
                );
            }
            else {
                str_ = to_string(name, derivative);
            }
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto str_state = to_string(
                SymEngine::down_cast<const TwoElecOperator&>(x).get_state()
            );
            auto derivative = SymEngine::down_cast<const TwoElecOperator&>(x).get_derivative();
            auto str_eri = derivative.empty()
                ? to_string(
                      std::string("ERI"),
                      SymEngine::down_cast<const TwoElecOperator&>(x).get_dependencies()
                  )
                : to_string(std::string("ERI"), derivative);
            auto name = SymEngine::down_cast<const TwoElecOperator&>(x).get_name();
            str_ = name + "(" + str_eri + ", " + str_state + ")";
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            str_ = to_string(
                SymEngine::down_cast<const ExchCorrPotential&>(x).get_name(),
                SymEngine::down_cast<const ExchCorrPotential&>(x).get_state(),
                SymEngine::down_cast<const ExchCorrPotential&>(x).get_derivative()
            );
        }
        else {
            throw SymEngine::NotImplementedError(
                "RemoveVisitor::bvisit() not implemented for " + x.__str__()
            );
        }
    }
}
