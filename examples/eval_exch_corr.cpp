#include <iostream>

#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"

using namespace Tinned;

void eval_xc_energy(const SymEngine::RCP<const ExchCorrEnergy>& energy)
{
    std::cout << "Output from eval_xc_energy()\n";
    auto state = energy->get_state();
    auto contr_map = energy->get_energy_map();
    for (const auto& term: contr_map) {
        for (const auto& contr: term.second) {
            std::cout << "Grid weight = " << stringify(term.first) << "\n";
            //get_derivative()
            std::cout << "XC functional order = " << contr.first << "\n";
            std::cout << "Unperturbed state = " << stringify(state) << "\n";
            if (!contr.second.is_null()) {
                std::cout << "Generalized density vector = "
                          << stringify(contr.second)
                          << "\n";
            }
        }
    }
    std::cout << "\n";
}

void eval_xc_potential(const SymEngine::RCP<const ExchCorrPotential>& potential)
{
    //auto = potential->get_potential_map();
    //auto vxc_map = extract_potential_map(potential);
}

// This is just a simple illustration, probably not for pratical use
int main()
{
    // Make perturbations
    auto a = make_perturbation(std::string("a"));
    auto b = make_perturbation(std::string("b"));
    auto c = make_perturbation(std::string("c"));
    auto d = make_perturbation(std::string("d"));

    // Make perturbations' dependecy and their maximum differentiated orders.
    // You can make different dependecies for different operators.
    auto dependencies = PertDependency({
        std::make_pair(a, 99),
        std::make_pair(b, 99),
        std::make_pair(c, 99),
        std::make_pair(d, 99)
    });

    // Make one-electron spin-orbital density matrix
    auto D = make_1el_density(std::string("D"));

    // Make grid weight
    auto weight = make_nonel_function(std::string("weight"), dependencies);

    // Make generalized overlap distribution
    auto Omega = make_1el_operator(std::string("Omega"), dependencies);

    // Make exchange-correlation energy functional and potential operator
    auto Exc = make_xc_energy(std::string("GGA"), D, Omega, weight);
    auto Vxc = make_xc_potential(std::string("GGA"), D, Omega, weight);

    // (0) Unperturbed case
    eval_xc_energy(Exc);
    eval_xc_potential(Vxc);

    // (1) The first order XC energy density derivative
    auto Exc_a = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc->diff(a));
    auto Vxc_a = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Vxc->diff(a));
    eval_xc_energy(Exc_a);
    eval_xc_potential(Vxc_a);

    // (2) The second order XC energy density derivative
    auto Exc_ab = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_a->diff(b));
    auto Vxc_ab = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Vxc_a->diff(b));
    eval_xc_energy(Exc_ab);
    eval_xc_potential(Vxc_ab);

    // (3) The third order XC energy density derivative
    auto Exc_abc = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_ab->diff(c));
    auto Vxc_abc = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Vxc_ab->diff(c));
    eval_xc_energy(Exc_abc);
    eval_xc_potential(Vxc_abc);

    // (4) The fourth order XC energy density derivative
    auto Exc_abcd = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(Exc_abc->diff(d));
    auto Vxc_abcd = SymEngine::rcp_dynamic_cast<const ExchCorrPotential>(Vxc_abc->diff(d));
    eval_xc_energy(Exc_abcd);
    eval_xc_potential(Vxc_abcd);

    return 0;
}
