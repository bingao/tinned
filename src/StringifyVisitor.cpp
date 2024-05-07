#include "Tinned/Perturbation.hpp"
#include "Tinned/LagMultiplier.hpp"

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

#include "Tinned/StateVector.hpp"
#include "Tinned/StateOperator.hpp"
#include "Tinned/AdjointMap.hpp"
#include "Tinned/ExpAdjointHamiltonian.hpp"

#include "Tinned/ZeroOperator.hpp"

#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    void StringifyVisitor::bvisit(const SymEngine::Symbol& x)
    {
        if (SymEngine::is_a_sub<const Perturbation>(x)) {
            auto& p = SymEngine::down_cast<const Perturbation&>(x);
            if (verbose_) {
                auto name = p.get_name();
                auto frequency = p.get_frequency();
                auto components = p.get_components();
                std::ostringstream o;
                o << apply(frequency) << "; ";
                if (components.empty()) {
                    o << "-";
                }
                else {
                    for (auto c=components.begin(); c!=components.end(); ++c) {
                        if (c==components.begin()) {
                            o << *c;
                        }
                        else {
                            o << "," << *c;
                        }
                    }
                }
                str_ = name + square_bracket(o.str());
            }
            else {
                str_ = p.get_name();
            }
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    void StringifyVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            auto& op = SymEngine::down_cast<const NonElecFunction&>(x);
            str_ = stringify_operator(
                op.get_name(), op.get_derivatives(), op.get_dependencies()
            );
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            auto str_eri = stringify_operator(
                std::string("ERI"), op.get_derivatives(), op.get_dependencies()
            );
            auto str_inner = stringify_state(op.get_inner_state());
            auto str_outer = stringify_state(op.get_outer_state());
            str_ = "tr" + square_bracket(
                "1/2*"+op.get_name()+parenthesize(str_eri+", "+str_inner)+"*"+str_outer
            );
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            auto order = op.get_order();
            auto str_op = order>0
                        ? op.get_name() + "^" + parenthesize(std::to_string(order))
                        : op.get_name();
            str_ = str_op + parenthesize(apply(op.get_inner()));
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            str_ = op.get_name() + parenthesize(apply(op.get_energy()));
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    void StringifyVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const LagMultiplier>(x)) {
            auto& op = SymEngine::down_cast<const LagMultiplier&>(x);
            str_ = stringify_operator(op.get_name(), op.get_derivatives());
        }
        else if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            auto& op = SymEngine::down_cast<const OneElecDensity&>(x);
            str_ = stringify_operator(op.get_name(), op.get_derivatives());
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const OneElecOperator&>(x);
            str_ = stringify_operator(
                op.get_name(), op.get_derivatives(), op.get_dependencies()
            );
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            auto str_eri = stringify_operator(
                std::string("ERI"), op.get_derivatives(), op.get_dependencies()
            );
            auto str_state = stringify_state(op.get_state());
            str_ = op.get_name() + parenthesize(str_eri+", "+str_state);
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            str_ = op.get_name() + parenthesize(apply(op.get_potential()));
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            str_ = op.get_name() + parenthesize(apply(op.get_target()));
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOverlap&>(x);
            str_ = op.get_name() + parenthesize(apply(op.get_braket()));
        }
        else if (SymEngine::is_a_sub<const ZeroOperator>(x)) {
            auto& op = SymEngine::down_cast<const ZeroOperator&>(x);
            str_ = stringify_operator(op.get_name());
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }
}
