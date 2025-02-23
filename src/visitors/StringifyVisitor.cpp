#include "Tinned/perturbations/Perturbation.hpp"

#include "Tinned/operators/PerturbedParameter.hpp"
#include "Tinned/operators/ConjugateTranspose.hpp"
#include "Tinned/operators/OneElecDensity.hpp"
#include "Tinned/operators/OneElecOperator.hpp"
#include "Tinned/operators/TwoElecEnergy.hpp"
#include "Tinned/operators/TwoElecOperator.hpp"
#include "Tinned/operators/CompositeFunction.hpp"
#include "Tinned/operators/ExchCorrEnergy.hpp"
#include "Tinned/operators/ExchCorrPotential.hpp"
#include "Tinned/operators/NonElecFunction.hpp"
#include "Tinned/operators/TemporumOperator.hpp"
#include "Tinned/operators/TemporumOverlap.hpp"
#include "Tinned/operators/AdjointMap.hpp"
#include "Tinned/operators/ClusterConjHamiltonian.hpp"
#include "Tinned/operators/ZeroOperator.hpp"

#include "Tinned/visitors/StringifyVisitor.hpp"

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
            std::ostringstream o;
            o << "1/2*" << apply(op.get_2el_operator())
              << "*" << stringify_state(op.get_outer_state());
            str_ = "tr" + square_bracket(o.str());
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            std::ostringstream o;
            o << op.get_name();
            auto order = op.get_order();
            if (order>0) o << "^" << parenthesize(std::to_string(order));
            str_ = o.str() + parenthesize(apply(op.get_inner()));
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            str_ = op.get_name() + parenthesize(apply(op.get_energy()));
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    void StringifyVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        if (SymEngine::is_a_sub<const ZeroOperator>(x)) {
            auto& op = SymEngine::down_cast<const ZeroOperator&>(x);
            str_ = stringify_operator(op.get_name());
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    void StringifyVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const PerturbedParameter>(x)) {
            auto& op = SymEngine::down_cast<const PerturbedParameter&>(x);
            str_ = stringify_operator(op.get_name(), op.get_derivatives());
        }
        else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            auto& op = SymEngine::down_cast<const ConjugateTranspose&>(x);
            str_ = op.get_name() + parenthesize(apply(op.get_arg()));
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
            std::ostringstream o;
            o << stringify_operator(std::string("ERI"), op.get_derivatives(), op.get_dependencies())
              << ", "
              << stringify_state(op.get_state());
            str_ = op.get_name() + parenthesize(o.str());
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
        else if (SymEngine::is_a_sub<const AdjointMap>(x)) {
            auto& op = SymEngine::down_cast<const AdjointMap&>(x);
            std::ostringstream o;
            for (const auto& term: op.get_x()) o << apply(term) + ", ";
            str_ = op.get_name() + parenthesize(o.str()+apply(op.get_y()));
        }
        else if (SymEngine::is_a_sub<const ClusterConjHamiltonian>(x)) {
            auto& op = SymEngine::down_cast<const ClusterConjHamiltonian&>(x);
            std::ostringstream o;
            o << apply(op.get_cluster_operator()) << ", " << apply(op.get_hamiltonian());
            str_ = op.get_name() + parenthesize(o.str());
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }
}
