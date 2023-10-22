#include <symengine/symengine_exception.h>

#include "Tinned/Perturbation.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/ElectronicState.hpp"
#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

#include "Tinned/ReplaceVisitor.hpp"

namespace Tinned
{
    void ReplaceVisitor::bvisit(const SymEngine::Symbol& x)
    {
        if (SymEngine::is_a_sub<const Perturbation>(x)) {
            auto name = SymEngine::down_cast<const Perturbation&>(x).get_name();
            auto dimension = SymEngine::down_cast<const Perturbation&>(x).get_dimension();
            std::ostringstream o;
            o << name << "[" << dimension << "]";
            str_ = o.str();
        }
        else {
            SymEngine::XReplaceVisitor::bvisit(x);
        }
    }

    void ReplaceVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            // We first
            for (const auto &p : subs_dict_) {
            }

            auto name = SymEngine::down_cast<const NonElecFunction&>(x).get_name();
            auto derivative = SymEngine::down_cast<const NonElecFunction&>(x).get_derivative();
            if (derivative.empty()) {
                str_ = to_string(
                    name,
                    SymEngine::down_cast<const NonElecFunction&>(x).get_dependencies()
                );
            }
            else {
                str_ = to_string(name, derivative);
            }
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            str_ = to_string(
                SymEngine::down_cast<const ExchCorrEnergy&>(x).get_name(),
                SymEngine::down_cast<const ExchCorrEnergy&>(x).get_state(),
                SymEngine::down_cast<const ExchCorrEnergy&>(x).get_derivative()
            );
        }
        else {
            SymEngine::XReplaceVisitor::bvisit(x);
        }
    }

    void bvisit(const SymEngine::ZeroMatrix& x);

    void ReplaceVisitor::bvisit(const SymEngine::MatrixSymbol& x)
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
            SymEngine::XReplaceVisitor::bvisit(x);
        }
    }

    void bvisit(const SymEngine::Trace& x);
    void bvisit(const SymEngine::ConjugateMatrix& x);
    void bvisit(const SymEngine::Transpose& x);
    void bvisit(const SymEngine::MatrixAdd& x);
    void bvisit(const SymEngine::MatrixMul& x);
    void bvisit(const SymEngine::MatrixDerivative& x);
}
