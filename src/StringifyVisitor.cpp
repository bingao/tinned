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
            auto derivative = op.get_derivative();
            if (derivative.empty()) {
                str_ = to_string(name, op.get_dependencies());
            }
            else {
                str_ = to_string(name, derivative);
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            str_ = to_string(op.get_name(), op.get_state(), op.get_derivative());
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
            auto derivative = op.get_derivative();
            str_ = derivative.empty() ? name : to_string(name, derivative);
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const OneElecOperator&>(x);
            auto name = op.get_name();
            auto derivative = op.get_derivative();
            str_ = derivative.empty()
                 ? to_string(name, op.get_dependencies())
                 : to_string(name, derivative);
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            auto str_state = to_string(op.get_state());
            auto derivative = op.get_derivative();
            auto str_eri = derivative.empty()
                ? to_string(std::string("ERI"), op.get_dependencies())
                : to_string(std::string("ERI"), derivative);
            str_ = op.get_name() + "(" + str_eri + ", " + str_state + ")";
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            str_ = to_string(op.get_name(), op.get_state(), op.get_derivative());
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
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    //std::string StringifyVisitor::apply(const SymEngine::Basic& x)
    //{
    //    x.accept(*this);
    //    return str_;
    //}

    std::string stringify(const SymEngine::Basic& x)
    {
        StringifyVisitor visitor;
        return visitor.apply(x);
    }
}
