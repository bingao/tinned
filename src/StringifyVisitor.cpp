#include "Tinned/Perturbation.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecEnergy.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

#include "Tinned/LagMultiplier.hpp"
#include "Tinned/StateVector.hpp"
#include "Tinned/StateOperator.hpp"
#include "Tinned/AdjointMap.hpp"
#include "Tinned/ExpAdjointHamiltonian.hpp"

#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    void StringifyVisitor::bvisit(const SymEngine::Symbol& x)
    {
        if (SymEngine::is_a_sub<const Perturbation>(x)) {
            auto& p = SymEngine::down_cast<const Perturbation&>(x);
            auto name = p.get_name();
            auto frequency = p.get_frequency();
            auto components = p.get_components();
            std::ostringstream o;
            o << name << "[" << apply(frequency) << "; ";
            if (components.empty()) {
                o << "-]";
            }
            else {
                for (auto c = components.begin(); c != components.end(); ++c) {
                    if (c == components.begin()) {
                        o << *c;
                    }
                    else {
                        o << "," << *c;
                    }
                }
                o << "]";
            }
            str_ = o.str();
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    void StringifyVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            auto& op = SymEngine::down_cast<const NonElecFunction&>(x);
            auto name = op.get_name();
            auto derivatives = op.get_derivatives();
            if (derivatives.empty()) {
                str_ = to_string(name, op.get_dependencies());
            }
            else {
                str_ = to_string(name, derivatives);
            }
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            auto str_inner = to_string(op.get_inner_state());
            auto str_outer = to_string(op.get_outer_state());
            auto derivatives = op.get_derivatives();
            auto str_eri = derivatives.empty()
                ? to_string(std::string("ERI"), op.get_dependencies())
                : to_string(std::string("ERI"), derivatives);
            str_ = "1/2*tr(" + op.get_name() + "(" + str_eri + ", " + str_inner
                 + ")*" + str_outer + ")";
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& fun = SymEngine::down_cast<const CompositeFunction&>(x);
            auto order = fun.get_order();
            auto str_fun = order > 0
                ? fun.get_name() + "^(" + std::to_string(order) + ")"
                : fun.get_name();
            str_ = str_fun + "(" + apply(*fun.get_inner()) + ")";
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            str_ = op.get_name() + "(" + apply(*op.get_energy()) + ")";
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    void StringifyVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            auto& op = SymEngine::down_cast<const OneElecDensity&>(x);
            auto name = op.get_name();
            auto derivatives = op.get_derivatives();
            str_ = derivatives.empty() ? name : to_string(name, derivatives);
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const OneElecOperator&>(x);
            auto name = op.get_name();
            auto derivatives = op.get_derivatives();
            str_ = derivatives.empty()
                 ? to_string(name, op.get_dependencies())
                 : to_string(name, derivatives);
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            auto str_state = to_string(op.get_state());
            auto derivatives = op.get_derivatives();
            auto str_eri = derivatives.empty()
                ? to_string(std::string("ERI"), op.get_dependencies())
                : to_string(std::string("ERI"), derivatives);
            str_ = op.get_name() + "(" + str_eri + ", " + str_state + ")";
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            str_ = op.get_name() + "(" + apply(*op.get_potential()) + ")";
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            std::ostringstream o;
            o << op.get_name() << "(" << apply(*op.get_target()) << ")";
            str_ = o.str();
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOverlap&>(x);
            std::ostringstream o;
            o << op.get_name() << "(" << apply(*op.get_braket()) << ")";
            str_ = o.str();
        }
        else if (SymEngine::is_a_sub<const LagMultiplier>(x)) {
            auto& op = SymEngine::down_cast<const LagMultiplier&>(x);
            auto name = op.get_name();
            auto derivatives = op.get_derivatives();
            str_ = derivatives.empty() ? name : to_string(name, derivatives);
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    //std::string StringifyVisitor::apply(const SymEngine::Basic& x)
    //{
    //    x.accept(*this);
    //    return str_;
    //}
}
