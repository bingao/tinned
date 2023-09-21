#include "Tinned/TinnedStrPrinter.hpp"

namespace Tinned
{
    void TinnedStrPrinter::bvisit(const SymEngine::Symbol& x)
    {
        if (SymEngine::is_a<const Perturbation>(x)) {
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

    void TinnedStrPrinter::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a<const OneElecOperator>(x)) {
            auto name = SymEngine::down_cast<const OneElecOperator&>(x).get_name();
            auto derivative = SymEngine::down_cast<const OneElecOperator&>(x).get_derivative();
            std::ostringstream o;
            if (derivative.empty()) {
                auto deps = SymEngine::down_cast<const OneElecOperator&>(x).get_dependencies();
                o << name << "(";
                for (auto p = deps.begin(); p != deps.end(); ++p) {
                    if (p == deps.begin()) {
                        o << apply(*p->first) << "<" << p->second << ">";
                    }
                    else {
                        o << ", " << apply(*p->first) << "<" << p->second << ">";
                    }
                }
                o << ")";
            }
            else {
                o << "Derivative(" << name;
                for (const auto& var: derivative) {
                    o << ", " << apply(*var);
                }
                o << ")";
            }
            str_ = o.str();
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
