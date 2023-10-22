#include <iostream>
#include <string>
#include <utility>

#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>

#include <symengine/subs.h>

#include <symengine/functions.h>
#include <symengine/mul.h>

#include "Tinned.hpp"

int main()
{
    auto f = SymEngine::make_rcp<const Tinned::Perturbation>(
        std::string("EL"),
        SymEngine::one
    );
    auto g = SymEngine::make_rcp<const Tinned::Perturbation>(
        std::string("GEO"),
        SymEngine::two
    );
    std::cout << Tinned::stringify(*f) << "\n";
    std::cout << Tinned::stringify(*g) << "\n\n";

    auto S = SymEngine::make_rcp<const Tinned::OneElecOperator>(
        std::string("S"),
        Tinned::PertDependency({
            //std::make_pair(f, 0),  //can be removed
            std::make_pair(g, 99)
        })
    );
    std::cout << Tinned::stringify(*S) << "\n\n";

    auto D = SymEngine::make_rcp<const Tinned::OneElecDensity>(std::string("D"));
    std::cout << Tinned::stringify(*D) << "\n\n";

    // idempotency constraint Z = DSD - D
    auto Z = SymEngine::matrix_add({
        SymEngine::matrix_mul({D, S, D}),
        SymEngine::matrix_mul({SymEngine::minus_one, D})
    });
    std::cout << Tinned::stringify(*Z) << "\n\n";

    auto K = (Z->diff(g))->diff(g);
    std::cout << Tinned::stringify(*K) << "\n\n";

    auto Dg = D->diff(g);
    auto Dgg = Dg->diff(g);
    auto Keq = Tinned::remove_if(K, SymEngine::vec_basic({Dg}));
    std::cout << "remove_if(): " << Tinned::stringify(*Keq) << "\n\n";
    Keq = Tinned::remove_if(K, SymEngine::vec_basic({Dgg}));
    std::cout << "remove_if(): " << Tinned::stringify(*Keq) << "\n\n";
    Keq = Tinned::remove_if(K, SymEngine::vec_basic({Dg, Dgg}));
    std::cout << "remove_if(): " << Tinned::stringify(*Keq) << "\n\n";

    auto Kneq = Tinned::keep_if(K, SymEngine::vec_basic({Dg}));
    std::cout << "keep_if(): " << Tinned::stringify(*Kneq) << "\n\n";
    Kneq = Tinned::keep_if(K, SymEngine::vec_basic({Dgg}));
    std::cout << "keep_if(): " << Tinned::stringify(*Kneq) << "\n\n";
    Kneq = Tinned::keep_if(K, SymEngine::vec_basic({Dg, Dgg}));
    std::cout << "keep_if(): " << Tinned::stringify(*Kneq) << "\n\n";

    // replacement
    //auto Dp = SymEngine::make_rcp<const Tinned::OneElecDensity>(std::string("Dp"));
    //SymEngine::map_basic_basic Msub;
    //Msub[D] = Dp;
    //auto M = Tinned::replace(K, Msub, true);
    //std::cout << "M = " << Tinned::stringify(*M) << "\n\n";

    auto r = SymEngine::function_symbol("r", SymEngine::vec_basic({f, g}));
    auto s = SymEngine::function_symbol("s", SymEngine::vec_basic({f, g}));
    auto prod = SymEngine::mul({r, s});
    //auto deriv = (prod->diff(f))->diff(g);
    auto deriv = prod->diff(f);
    SymEngine::map_basic_basic psub;
    psub[r] = f;
    psub[g] = f;
    //psub[(r->diff(f))->diff(f)] = f;
    std::cout << "prod = " << Tinned::stringify(*prod) << "\n";
    std::cout << "deriv = " << Tinned::stringify(*deriv) << "\n";
    std::cout << "msubs = " << Tinned::stringify(*SymEngine::msubs(deriv, psub)) << "\n";
    std::cout << "xreplace = " << Tinned::stringify(*deriv->xreplace(psub)) << "\n\n";

    return 0;

    auto G = SymEngine::make_rcp<const Tinned::TwoElecOperator>(
        std::string("G"),
        D,
        Tinned::PertDependency({std::make_pair(g, 99)})
    );
    std::cout << Tinned::stringify(*G) << "\n\n";
    std::cout << Tinned::stringify(*(((G->diff(g))->diff(f))->diff(g))) << "\n\n";

    auto Exc = SymEngine::make_rcp<const Tinned::ExchCorrEnergy>(
        std::string("Exc"), D
    );
    std::cout << Tinned::stringify(*Exc) << "\n\n";
    std::cout << Tinned::stringify(*(((Exc->diff(g))->diff(f))->diff(g))) << "\n\n";

    auto Fxc = SymEngine::make_rcp<const Tinned::ExchCorrPotential>(
        std::string("Fxc"), D
    );
    std::cout << Tinned::stringify(*Fxc) << "\n\n";
    std::cout << Tinned::stringify(*(((Fxc->diff(g))->diff(f))->diff(g))) << "\n\n";

    auto h_nuc = SymEngine::make_rcp<const Tinned::NonElecFunction>(
        std::string("h_nuc"),
        Tinned::PertDependency({
            std::make_pair(f, 99),
            std::make_pair(g, 99)
        })
    );
    std::cout << Tinned::stringify(*h_nuc) << "\n\n";
    std::cout << Tinned::stringify(*(((h_nuc->diff(g))->diff(f))->diff(g))) << "\n\n";

    auto Dt = SymEngine::make_rcp<const Tinned::TemporumOperator>(D);
    auto Dw = SymEngine::rcp_dynamic_cast<const Tinned::TemporumOperator>((Dt->diff(g))->diff(g));
    std::cout << Tinned::stringify(*Dt) << "\n\n";
    std::cout << Tinned::stringify(*Dw) << "\n";
    std::cout << "frequnecy: " << Tinned::stringify(*Dw->get_frequency()) << "\n\n";

    auto T = SymEngine::make_rcp<const Tinned::TemporumOverlap>(
        Tinned::PertDependency({std::make_pair(g, 99), std::make_pair(f, 99)})
    );
    std::cout << Tinned::stringify(*T) << "\n\n";
    auto Tg = T->diff(g);
    std::cout << Tinned::stringify(*Tg) << "\n\n";

    auto Tgg = SymEngine::rcp_dynamic_cast<const Tinned::TemporumOverlap>(Tg->diff(g));
    std::cout << Tinned::stringify(*Tgg) << "\n";
    std::cout << "size: " << Tgg->size() << "\n";
    for (std::size_t i = 0; i<Tgg->size(); ++i) {
        std::cout << "frequnecy: " << Tinned::stringify(*Tgg->get_frequency(i));
        auto braket = Tgg->get_braket_product(i);
        std::cout << ", coef.: " << Tinned::stringify(*std::get<0>(braket))
                  << ", bra: " << Tinned::stringify(*std::get<1>(braket))
                  << ", ket: " << Tinned::stringify(*std::get<2>(braket))
                  << "\n\n";
    }
    auto braket1 = Tgg->get_braket_product(1);
    auto braket2 = Tgg->get_braket_product(2);
    std::cout << "equal: " << std::get<0>(braket1)->__eq__(*std::get<0>(braket2)) << "\n";
    std::cout << "equal: " << std::get<1>(braket1)->__eq__(*std::get<1>(braket2)) << "\n";
    std::cout << "bra: " << Tinned::stringify(*std::get<1>(braket1)) << "\n"
              << "bra: " << Tinned::stringify(*std::get<1>(braket2)) << "\n";
    std::cout << "equal: " << std::get<2>(braket1)->__eq__(*std::get<2>(braket2)) << "\n";
    std::cout << "ket: " << Tinned::stringify(*std::get<2>(braket1)) << "\n"
              << "ket: " << Tinned::stringify(*std::get<2>(braket2)) << "\n";

    return 0;
}
