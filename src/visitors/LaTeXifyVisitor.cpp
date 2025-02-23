#include <symengine/constants.h>

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

#include "Tinned/visitors/LaTeXifyVisitor.hpp"

namespace Tinned
{
    void LaTeXifyVisitor::bvisit(const SymEngine::Symbol& x)
    {
        if (SymEngine::is_a_sub<const Perturbation>(x)) {
            auto& p = SymEngine::down_cast<const Perturbation&>(x);
            str_ = p.get_name();
        }
        else {
            SymEngine::LatexPrinter::bvisit(x);
        }
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::FunctionSymbol& x)
    {
        if (SymEngine::is_a_sub<const NonElecFunction>(x)) {
            auto& op = SymEngine::down_cast<const NonElecFunction&>(x);
            str_ = latexify_operator(
                op.get_name(), op.get_derivatives(), OperFontStyle::Regular
            );
            update_num_symbols(1, str_);
        }
        else if (SymEngine::is_a_sub<const TwoElecEnergy>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecEnergy&>(x);
            std::ostringstream o;
            // `print_mul` is necessary, in particular if `apply` ends with a newline
            o << "\\tfrac{1}{2}"
              << apply(op.get_2el_operator())
              << print_mul()
              << apply(op.get_outer_state());
            str_ = "\\mathrm{tr}" + parenthesize(o.str());
        }
        else if (SymEngine::is_a_sub<const CompositeFunction>(x)) {
            auto& op = SymEngine::down_cast<const CompositeFunction&>(x);
            std::ostringstream o;
            o << op.get_name();
            auto order = op.get_order();
            if (order>0) o << "^{" << parenthesize(std::to_string(order)) << "}";
            str_ = o.str() + parenthesize(apply(op.get_inner()));
        }
        else if (SymEngine::is_a_sub<const ExchCorrEnergy>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrEnergy&>(x);
            str_ = latexify_operator(
                op.get_name(), op.get_derivatives(), OperFontStyle::Regular
            );
            update_num_symbols(1, str_);
        }
        else {
            SymEngine::LatexPrinter::bvisit(x);
            //FIXME: only one symbol?
            update_num_symbols(1, str_);
        }
    }

    // This function is modified from that of `symengine/printers/strprinter.cpp`
    void LaTeXifyVisitor::bvisit(const SymEngine::Add& x)
    {
        std::ostringstream o;
        bool first = true;
        SymEngine::map_basic_num dict(x.get_dict().begin(), x.get_dict().end());
        // Parentheses will be needed for more than one terms addition
        std::size_t num_terms = 0;
        if (SymEngine::neq(*(x.get_coef()), *SymEngine::zero)) {
            o << apply(x.get_coef());
            first = false;
            ++num_terms;
        }
        for (const auto& p: dict) {
            zero_factor_ = false;
            auto curr_num_symbols = curr_num_symbols_;
            std::string t;
            if (SymEngine::eq(*(p.second), *SymEngine::one)) {
                t = parenthesizeLT(p.first, SymEngine::PrecedenceEnum::Add);
            }
            else if (SymEngine::eq(*(p.second), *SymEngine::minus_one)) {
                t = "-" + parenthesizeLT(p.first, SymEngine::PrecedenceEnum::Mul);
            }
            else {
                // `p.second` cannot be zero, so we latexify `p.first` at the end
                t = parenthesizeLT(p.second, SymEngine::PrecedenceEnum::Mul)
                  + print_mul();
                t += parenthesizeLT(p.first, SymEngine::PrecedenceEnum::Mul);
            }
            if (zero_factor_) {
                curr_num_symbols_ = curr_num_symbols;
                continue;
            }
            ++num_terms;
            if (not first) {
                if (t[0]=='-') {
                    // `Trace` may have -1 as the coefficient factor
                    if (t.substr(0, 3)==std::string("-1 ")) {
                        o << "-" << t.substr(3);
                    }
                    else {
                        o << t;
                    }
                }
                else {
                    o << "+" << t;
                }
            }
            else {
                o << t;
                first = false;
            }
        }
        if (num_terms==0) {
            str_ = "";
        }
        else {
            // Important to set `zero_factor_` in particular the last term is zero
            zero_factor_ = false;
            str_ = o.str();
            //if (num_terms==1) {
            //    str_ = o.str();
            //}
            //else {
            //    str_ = parenthesize(o.str());
            //}
        }
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::Mul& x)
    {
        zero_factor_ = false;
        auto curr_num_symbols = curr_num_symbols_;
        SymEngine::LatexPrinter::bvisit(x);
        if (zero_factor_) {
            str_.clear();
            curr_num_symbols_ = curr_num_symbols;
        }
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::ZeroMatrix& x)
    {
        if (SymEngine::is_a_sub<const ZeroOperator>(x)) {
            auto& op = SymEngine::down_cast<const ZeroOperator&>(x);
            str_ = latexify_operator(op.get_name());
        }
        else {
            SymEngine::LatexPrinter::bvisit(x);
        }
        update_num_symbols(1, str_);
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::MatrixSymbol& x)
    {
        if (SymEngine::is_a_sub<const PerturbedParameter>(x)) {
            auto& op = SymEngine::down_cast<const PerturbedParameter&>(x);
            str_ = latexify_operator(op.get_name(), op.get_derivatives());
            update_num_symbols(1, str_);
        }
        else if (SymEngine::is_a_sub<const ConjugateTranspose>(x)) {
            auto op = SymEngine::down_cast<const ConjugateTranspose&>(x).get_arg();
            str_ = SymEngine::is_a_sub<const SymEngine::MatrixSymbol>(*op)
                 ? add_suffix(apply(op), "^{\\dagger}")
                 : add_suffix(parenthesize(apply(op)), "^{\\dagger}");
        }
        else if (SymEngine::is_a_sub<const OneElecDensity>(x)) {
            auto& op = SymEngine::down_cast<const OneElecDensity&>(x);
            str_ = latexify_operator(op.get_name(), op.get_derivatives());
            update_num_symbols(1, str_);
        }
        else if (SymEngine::is_a_sub<const OneElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const OneElecOperator&>(x);
            str_ = latexify_operator(op.get_name(), op.get_derivatives());
            update_num_symbols(1, str_);
        }
        else if (SymEngine::is_a_sub<const TwoElecOperator>(x)) {
            auto& op = SymEngine::down_cast<const TwoElecOperator&>(x);
            auto str_op = latexify_operator(op.get_name(), op.get_derivatives()) + "(";
            update_num_symbols(1, str_op);
            str_ = str_op + add_suffix(apply(op.get_state()), ")");
        }
        else if (SymEngine::is_a_sub<const ExchCorrPotential>(x)) {
            auto& op = SymEngine::down_cast<const ExchCorrPotential&>(x);
            str_ = latexify_operator(op.get_name(), op.get_derivatives());
            update_num_symbols(1, str_);
        }
        else if (SymEngine::is_a_sub<const TemporumOperator>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOperator&>(x);
            std::ostringstream o;
            if (op.get_type()==TemporumType::Bra) o << "-";
            o << "\\text{i}\\frac{\\partial}{\\partial t}" << apply(op.get_target());
            str_ = o.str();
        }
        else if (SymEngine::is_a_sub<const TemporumOverlap>(x)) {
            auto& op = SymEngine::down_cast<const TemporumOverlap&>(x);
            auto derivatives = op.get_derivatives();
            if (derivatives.empty()) {
                // Unperturbed T matrix is actually a zero operator
                str_ = latexify_operator(op.get_name()) + "^{0}";
                update_num_symbols(1, str_);
            }
            else {
                str_ = latexify_operator(op.get_name(), derivatives);
                update_num_symbols(1, str_);
            }
        }
        else if (SymEngine::is_a_sub<const AdjointMap>(x)) {
            auto& op = SymEngine::down_cast<const AdjointMap&>(x);
            LaTeXifyVisitor visitor;
            auto terms = op.get_x();
            std::string str_x;
            for (std::size_t i=0; i<terms.size(); ++i) {
                if (i>0) str_x += print_mul();
                str_x += "(\\mathrm{" + op.get_name() + "}_{"
                       + remove_newline(visitor.apply(terms[i])) + "})";
                if (i==terms.size()-1) str_x += "(";
                update_num_symbols(1, str_x);
            }
            str_ = str_x + add_suffix(apply(op.get_y()), ")");
        }
        else if (SymEngine::is_a_sub<const ClusterConjHamiltonian>(x)) {
            auto& op = SymEngine::down_cast<const ClusterConjHamiltonian&>(x);
            auto cluster_op = op.get_cluster_operator();
            // Subscripts and superscripts should not change the number of
            // symbols, so we use a new visitor
            LaTeXifyVisitor visitor;
            auto str_cluster_op
                = SymEngine::is_a_sub<const SymEngine::MatrixSymbol>(*cluster_op)
                ? visitor.apply(cluster_op) : parenthesize(visitor.apply(cluster_op));
            auto str_op = "\\mathrm{e}^{\\mathrm{ad}_{-"
                        + remove_newline(str_cluster_op) + "}}(";
            update_num_symbols(2, str_op);
            str_ = str_op + add_suffix(apply(op.get_hamiltonian()), ")");
        }
        else {
            SymEngine::LatexPrinter::bvisit(x);
            //FIXME: only one symbol?
            update_num_symbols(1, str_);
        }
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::MatrixAdd& x)
    {
        std::ostringstream o;
        auto terms = x.get_terms();
        auto num_terms = terms.size();
        for (std::size_t i=0; i<terms.size(); ++i) {
            zero_factor_ = false;
            auto curr_num_symbols = curr_num_symbols_;
            auto str_term = apply(terms[i]);
            if (zero_factor_) {
                curr_num_symbols_ = curr_num_symbols;
                --num_terms;
            }
            else {
                if (i>0) o << "+";
                o << str_term;
            }
        }
        if (num_terms==0) {
            str_ = "";
        }
        else {
            // Important to set `zero_factor_` in particular the last term is zero
            zero_factor_ = false;
            if (num_terms==1) {
                str_ = o.str();
                // `MatrixMul` may have -1 as the coefficient factor
                if (str_.substr(0, 3)==std::string("-1 ")) str_.erase(1, 2);
            }
            else {
                str_ = parenthesize(o.str());
                str_ = std::regex_replace(str_, std::regex(R"(-1\s)"), R"(-)");
                str_ = std::regex_replace(str_, std::regex(R"(\+-)"), R"(-)");
            }
        }
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::MatrixMul& x)
    {
        zero_factor_ = false;
        auto curr_num_symbols = curr_num_symbols_;
        SymEngine::LatexPrinter::bvisit(x);
        if (zero_factor_) {
            str_.clear();
            curr_num_symbols_ = curr_num_symbols;
        }
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::Trace& x)
    {
        str_ = "\\mathrm{tr}" + parenthesize(apply(x.get_args()[0]));
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::ConjugateMatrix& x)
    {
        auto op = x.get_arg();
        str_ = SymEngine::is_a_sub<const SymEngine::MatrixSymbol>(*op)
             ? add_suffix(apply(op), "^{*}")
             : add_suffix(parenthesize(apply(op)), "^{*}");
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::Transpose& x)
    {
        auto op = x.get_arg();
        str_ = SymEngine::is_a_sub<const SymEngine::MatrixSymbol>(*op)
             ? add_suffix(apply(op), "^{\\mathrm{T}}")
             : add_suffix(parenthesize(apply(op)), "^{\\mathrm{T}}");
    }

    void LaTeXifyVisitor::bvisit(const SymEngine::MatrixDerivative& x)
    {
        auto op = x.get_arg();
        str_ = SymEngine::is_a_sub<const SymEngine::MatrixSymbol>(*op)
             ? add_suffix(apply(op), latexify_derivatives(x.get_symbols()))
             : add_suffix(parenthesize(apply(op)), latexify_derivatives(x.get_symbols()));
    }

    std::string LaTeXifyVisitor::parenthesize(const std::string& expr)
    {
        return "(" + add_suffix(expr, ")");
    }

    std::string LaTeXifyVisitor::print_mul()
    {
        return find_trailing_newline(str_)==std::string::npos ? " " : "\\times ";
    }
}
