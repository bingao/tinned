#include <iostream>

#include "Tinned/TinnedStrPrinter.hpp"

namespace Tinned
{
    void TinnedStrPrinter::bvisit(const SymEngine::Symbol& x)
    {
        if (SymEngine::is_a_sub<const Perturbation>(x)) {
            auto name = SymEngine::down_cast<const Perturbation&>(x).get_name();
            auto dimension = SymEngine::down_cast<const Perturbation&>(x).get_dimension();
            std::ostringstream o;
            o << name << "[" << dimension << "]";
            str_ = o.str();
        }
        else {
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    void TinnedStrPrinter::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
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
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    void TinnedStrPrinter::bvisit(const SymEngine::MatrixSymbol& x)
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
            SymEngine::StrPrinter::bvisit(x);
        }
    }

    //std::string TinnedStrPrinter::apply(const SymEngine::Basic& x)
    //{
    //    x.accept(*this);
    //    return str_;
    //}

    std::string tinned_str(const SymEngine::Basic& x)
    {
        TinnedStrPrinter strPrinter;
        return strPrinter.apply(x);
    }
}
