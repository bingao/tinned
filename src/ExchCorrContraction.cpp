#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/number.h>
#include <symengine/constants.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrContraction.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/KeepVisitor.hpp"
#include "Tinned/StringifyVisitor.hpp"

namespace Tinned
{
    std::tuple<SymEngine::RCP<const NonElecFunction>,
               SymEngine::RCP<const CompositeFunction>,
               SymEngine::RCP<const SymEngine::Basic>>
    extract_exc_contraction(const SymEngine::RCP<const SymEngine::Mul>& expression)
    {
        SymEngine::RCP<const NonElecFunction> weight;
        SymEngine::RCP<const CompositeFunction> exc;
        // Hold generalized density vectors and/or coefficient
        SymEngine::vec_basic factors;
        for (const auto& arg: expression->get_args()) {
            // Coefficient
            if (SymEngine::is_a_sub<const SymEngine::Number>(*arg)) {
                factors.push_back(arg);
            }
            // Grid weight
            else if (SymEngine::is_a_sub<const NonElecFunction>(*arg)) {
                if (!weight.is_null()) throw SymEngine::SymEngineException(
                    "Two grid weights got from the XC energy density contraction "
                    + stringify(expression)
                );
                weight = SymEngine::rcp_dynamic_cast<const NonElecFunction>(arg);
            }
            // XC functional derivative
            else if (SymEngine::is_a_sub<const CompositeFunction>(*arg)) {
                if (!exc.is_null()) throw SymEngine::SymEngineException(
                    "Two XC functional derivatives got from the XC energy density contraction "
                    + stringify(expression)
                );
                exc = SymEngine::rcp_dynamic_cast<const CompositeFunction>(arg);
            }
            // Generalized density vector, sum of generalized density
            // vectors, or power of (sum of) generalized density
            // vector(s)
            else if (
                SymEngine::is_a_sub<const SymEngine::Trace>(*arg) ||
                SymEngine::is_a_sub<const SymEngine::Add>(*arg) ||
                SymEngine::is_a_sub<const SymEngine::Pow>(*arg)
            ) {
                factors.push_back(arg);
            }
            else {
                throw SymEngine::SymEngineException(
                    "Invalid type from the XC energy density contraction "
                    + stringify(expression)
                );
            }
        }
        if (weight.is_null() || exc.is_null()) {
            throw SymEngine::SymEngineException(
                "Terms missing ("
                + std::to_string(weight.is_null()) + "/"
                + std::to_string(exc.is_null()) +
                + ") from the XC energy density contraction "
                + stringify(expression)
            );
        }
        else {
            if (factors.empty()) {
                // We return nullptr when there is no generalized density
                // vectors, for example in the unperturbed case
                return std::make_tuple(
                    weight, exc, SymEngine::RCP<const SymEngine::Basic>()
                );
            }
            else {
                auto dens_vectors = SymEngine::mul(factors);
                // We need to make sure there is neither grid weights nor XC
                // functional derivatives in the generalized density vectors
                auto factors_left = keep_if(
                    dens_vectors, SymEngine::set_basic({weight, exc})
                );
                if (factors_left.is_null()) {
                    return std::make_tuple(weight, exc, dens_vectors);
                }
                else {
                    throw SymEngine::SymEngineException(
                        "Invalid grid weights and/or XC functional derivatives "
                        + stringify(factors_left)
                        + " in "
                        + stringify(dens_vectors)
                    );
                }
            }
        }
    }

    void add_exc_contraction(
        ExcContractionMap& energyMap,
        const std::tuple<SymEngine::RCP<const NonElecFunction>,
                         SymEngine::RCP<const CompositeFunction>,
                         SymEngine::RCP<const SymEngine::Basic>>& value
    )
    {
        // Check if the grid weight exists in the outer map
        auto weight_map = energyMap.find(std::get<0>(value));
        if (weight_map==energyMap.end()) {
            energyMap.emplace(
                std::get<0>(value),
                ExcDensityContractionMap({{std::get<1>(value), std::get<2>(value)}})
            );
        }
        else {
            // Check if the functional derivative of XC energy density exists
            // in the inner map
            auto exc_map = weight_map->second.find(std::get<1>(value));
            if (exc_map==weight_map->second.end()) {
                weight_map->second.emplace(std::get<1>(value), std::get<2>(value));
            }
            else {
                exc_map->second = SymEngine::add(exc_map->second, std::get<2>(value));
            }
        }
    }

    // Extract all terms in a sum of XC energy or its derivatives
    ExcContractionMap extract_energy_map_(
        const SymEngine::RCP<const SymEngine::Add>& expression,
        const SymEngine::RCP<const SymEngine::Number>& coef = SymEngine::one
    )
    {
        ExcContractionMap energy_map;
        auto sum_energies = expression->get_args();
        // No coefficient exists in the XC energy derivatives
        SYMENGINE_ASSERT(
            !SymEngine::is_a_sub<const SymEngine::Number>(*sum_energies.front())
        )
        for (const auto& expr_energy: sum_energies) {
            SYMENGINE_ASSERT(
                SymEngine::is_a_sub<const SymEngine::Mul>(*expr_energy)
            )
            auto energy = SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(expr_energy);
            auto args = energy->get_args();
            // This is another sum of XC energy or its derivatives
            if (args.size()==2
                && SymEngine::is_a_sub<const SymEngine::Number>(*args.front())
                && SymEngine::is_a_sub<const SymEngine::Add>(*args.back())) {
                if (energy_map.empty()) {
                    energy_map = extract_energy_map_(
                        SymEngine::rcp_dynamic_cast<const SymEngine::Add>(args.back()),
                        SymEngine::rcp_dynamic_cast<const SymEngine::Number>(args.front())
                    );
                }
                else {
                    merge_energy_map(
                        energy_map,
                        extract_energy_map_(
                            SymEngine::rcp_dynamic_cast<const SymEngine::Add>(args.back()),
                            SymEngine::rcp_dynamic_cast<const SymEngine::Number>(args.front())
                        )
                    );
                }
            }
            else {
                auto contr_term = coef->is_one()
                    ? extract_exc_contraction(
                          SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(expr_energy)
                      )
                    : extract_exc_contraction(
                          SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(
                              SymEngine::mul(coef, expr_energy)
                          )
                      );
                add_exc_contraction(energy_map, contr_term);
            }
        }
        return energy_map;
    }

    ExcContractionMap extract_energy_map(
        const SymEngine::RCP<const SymEngine::Basic>& expression
    )
    {
        // We first expand `expression` so that generalized density vectors
        // could be collected and simplified
        auto expr_expand = SymEngine::expand(expression);
        // XC energy or its derivatives must be either `SymEngine::Mul` or
        // `SymEngine::Add`
        SYMENGINE_ASSERT(
            SymEngine::is_a<const SymEngine::Mul>(*expr_expand) ||
            SymEngine::is_a<const SymEngine::Add>(*expr_expand)
        )
        // Unperturbed, the first-order case or `SymEngine::Add` multiplied by
        // a coefficient (e.g. results from `keep_if()`)
        if (SymEngine::is_a_sub<const SymEngine::Mul>(*expr_expand)) {
            auto energy = SymEngine::rcp_dynamic_cast<const SymEngine::Mul>(expr_expand);
            auto args = energy->get_args();
            if (args.size()==2
                && SymEngine::is_a_sub<const SymEngine::Number>(*args.front())
                && SymEngine::is_a_sub<const SymEngine::Add>(*args.back()))
            {
                return extract_energy_map_(
                    SymEngine::rcp_dynamic_cast<const SymEngine::Add>(args.back()),
                    SymEngine::rcp_dynamic_cast<const SymEngine::Number>(args.front())
                );
            }
            else {
                auto contr_term = extract_exc_contraction(energy);
                return ExcContractionMap({
                    {
                        std::get<0>(contr_term),
                        ExcDensityContractionMap({
                            {std::get<1>(contr_term), std::get<2>(contr_term)}
                        })
                    }
                });
            }
        }
        // Perturbed case
        else {
            return extract_energy_map_(
                SymEngine::rcp_dynamic_cast<const SymEngine::Add>(expr_expand)
            );
        }
    }

    void merge_energy_map(ExcContractionMap& map1, const ExcContractionMap& map2)
    {
        for (const auto& weight_map2: map2) {
            auto weight_map1 = map1.find(weight_map2.first);
            // New grid weight and corresponding contractions
            if (weight_map1==map1.end()) {
                 map1.emplace(weight_map2);
            }
            else {
                for (const auto& exc_map2: weight_map2.second) {
                    auto exc_map1 = weight_map1->second.find(exc_map2.first);
                    // New XC energy derivative found
                    if (exc_map1==weight_map1->second.end()) {
                        weight_map1->second.emplace(exc_map2);
                    }
                    else {
                        exc_map1->second = SymEngine::add(
                            exc_map1->second, exc_map2.second
                        );
                    }
                }
            }
        }
    }

    SymEngine::RCP<const SymEngine::Basic> convert_energy_map(
        const ExcContractionMap& energyMap
    )
    {
        SymEngine::vec_basic terms = {};
        // `weight_map.first` is the (un)perturbed grid weight, and
        // `weight_map.second` is type `ExcDensityContractionMap`
        for (const auto& weight_map: energyMap) {
            // `exc_map.first` is the order of XC energy functional derivative
            // vector(s), and `exc_map.second` is the generalized density
            // vector(s)
            for (const auto& exc_map: weight_map.second) {
                if (exc_map.second.is_null()) {
                    terms.push_back(SymEngine::mul(weight_map.first, exc_map.first));
                }
                else {
                    terms.push_back(SymEngine::mul(SymEngine::vec_basic({
                        weight_map.first, exc_map.first, exc_map.second
                    })));
                }
            }
        }
        SYMENGINE_ASSERT(!terms.empty())
        return terms.size()==1 ? terms[0] : SymEngine::add(terms);
    }

    bool eq_energy_map(
        const ExcContractionMap& map1, const ExcContractionMap& map2
    )
    {
        if (map1.size()==map2.size()) {
            auto weight_map1 = map1.begin();
            auto weight_map2 = map2.begin();
            for (; weight_map1!=map1.end(); ++weight_map1, ++weight_map2) {
                // Compare grid weights
                if (SymEngine::eq(*weight_map1->first, *weight_map2->first)) {
                    if (weight_map1->second.size()==weight_map2->second.size()) {
                        auto exc_map1 = weight_map1->second.begin();
                        auto exc_map2 = weight_map2->second.begin();
                        for (; exc_map1!=weight_map1->second.end();) {
                            if (SymEngine::eq(*exc_map1->first, *exc_map2->first)) {
                                if (exc_map1->second.is_null()) {
                                    if (!exc_map2->second.is_null()) return false;
                                }
                                else {
                                    if (exc_map2->second.is_null()) return false;
                                    if (SymEngine::neq(
                                        *exc_map1->second, *exc_map2->second
                                    )) return false;
                                }
                            }
                            else {
                                return false;
                            }
                            ++exc_map1;
                            ++exc_map2;
                        }
                    }
                    else {
                        return false;
                    }
                }
                else {
                    return false;
                }
            }
            return true;
        }
        else {
            return false;
        }
    }

    std::pair<SymEngine::RCP<const OneElecOperator>, ExcContractionMap>
    extract_vxc_contraction(
        const SymEngine::RCP<const SymEngine::MatrixMul>& expression
    )
    {
        auto energy_map = extract_energy_map(expression->get_scalar());
        auto factors = expression->get_factors();
        if (factors.size()==1) {
            SYMENGINE_ASSERT(SymEngine::is_a_sub<const OneElecOperator>(*factors[0]))
            auto Omega = SymEngine::rcp_dynamic_cast<const OneElecOperator>(factors[0]);
            if (energy_map.empty()) {
                throw SymEngine::SymEngineException(
                    "Empty ExcContractionMap from " + stringify(expression)
                );
            }
            else {
                return std::make_pair(Omega, energy_map);
            }
        }
        else {
            throw SymEngine::SymEngineException(
                "Invalid overlap distributions in the XC potential contraction "
                + stringify(expression)
            );
        }
    }

    VxcContractionMap extract_potential_map(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& expression
    )
    {
        // XC potential or its derivatives must be either
        // `SymEngine::MatrixMul` or `SymEngine::MatrixAdd`
        SYMENGINE_ASSERT(
            SymEngine::is_a<const SymEngine::MatrixMul>(*expression) ||
            SymEngine::is_a<const SymEngine::MatrixAdd>(*expression)
        )
        // Unperturbed case or when the generalized overlap distribution does
        // not depend on the applied perturbation(s)
        if (SymEngine::is_a_sub<const SymEngine::MatrixMul>(*expression)) {
            auto potential_term = extract_vxc_contraction(
                SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(expression)
            );
            return VxcContractionMap({potential_term});
        }
        // Perturbed case
        else {
            VxcContractionMap potential_map;
            auto potential = SymEngine::rcp_dynamic_cast<const SymEngine::MatrixAdd>(expression);
            auto potential_terms = potential->get_args();
            for (const auto& term: potential_terms) {
                SYMENGINE_ASSERT(
                    SymEngine::is_a_sub<const SymEngine::MatrixMul>(*term)
                )
                auto vxc_factors = extract_vxc_contraction(
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixMul>(term)
                );
                // Check if the generalized overlap distribution exists
                // in the outermost map
                auto energy_map = potential_map.find(vxc_factors.first);
                if (energy_map==potential_map.end()) {
                    potential_map.emplace(vxc_factors);
                }
                else {
                    merge_energy_map(energy_map->second, vxc_factors.second);
                }
            }
            return potential_map;
        }
    }

    void merge_potential_map(VxcContractionMap& map1, const VxcContractionMap& map2)
    {
        for (const auto& energy_map2: map2) {
            auto energy_map1 = map1.find(energy_map2.first);
            // New generalized overlap distribution
            if (energy_map1==map1.end()) {
                 map1.emplace(energy_map2);
            }
            else {
                merge_energy_map(energy_map1->second, energy_map2.second);
            }
        }
    }

    SymEngine::RCP<const SymEngine::MatrixExpr> convert_potential_map(
        const VxcContractionMap& potentialMap
    )
    {
        SymEngine::vec_basic factors = {};
        // `energy_map.first` is the (un)perturbed generalized overlap
        // distribution, and `energy_map.second` is type `ExcContractionMap`
        for (const auto& energy_map: potentialMap) {
            factors.push_back(
                SymEngine::matrix_mul(SymEngine::vec_basic({
                    convert_energy_map(energy_map.second),
                    energy_map.first
                }))
            );
        }
        SYMENGINE_ASSERT(!factors.empty())
        return factors.size()==1
            ? SymEngine::rcp_dynamic_cast<const SymEngine::MatrixExpr>(factors[0])
            : SymEngine::matrix_add(factors);
    }

    bool eq_potential_map(const VxcContractionMap& map1, const VxcContractionMap& map2)
    {
        if (map1.size()==map2.size()) {
            auto energy_map1 = map1.begin();
            auto energy_map2 = map2.begin();
            for (; energy_map1!=map1.end(); ++energy_map1, ++energy_map2) {
                // Compare generalized overlap distribution
                if (SymEngine::eq(*energy_map1->first, *energy_map2->first)) {
                    if (!eq_energy_map(energy_map1->second, energy_map2->second))
                        return false;
                }
                else {
                    return false;
                }
            }
            return true;
        }
        else {
            return false;
        }
    }
}
