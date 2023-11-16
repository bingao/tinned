#include <string>

#include <symengine/pow.h>
#include <symengine/number.h>
#include <symengine/constants.h>
#include <symengine/matrices/trace.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrContraction.hpp"
#include "Tinned/CompositeFunction.hpp"
#include "Tinned/KeepVisitor.hpp"

namespace Tinned
{
    std::tuple<SymEngine::RCP<const NonElecFunction>,
               unsigned int,
               SymEngine::RCP<const SymEngine::Basic>>
    extract_exc_contraction(const SymEngine::RCP<const SymEngine::Mul>& expression)
    {
        SymEngine::RCP<const NonElecFunction> weight;
        SymEngine::RCP<const CompositeFunction> exc;
        unsigned int order;
        bool found_order = false;
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
                    + expression->__str__()
                );
                weight = SymEngine::rcp_dynamic_cast<const NonElecFunction>(arg);
            }
            // XC functional derivative
            else if (SymEngine::is_a_sub<const CompositeFunction>(*arg)) {
                if (found_order) throw SymEngine::SymEngineException(
                    "Two XC functional derivatives got from the XC energy density contraction "
                    + expression->__str__()
                );
                exc = SymEngine::rcp_dynamic_cast<const CompositeFunction>(arg);
                order = exc->get_order();
                found_order = true;
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
                    + expression->__str__()
                );
            }
        }
        if (weight.is_null() || !found_order) {
            throw SymEngine::SymEngineException(
                "Terms missing ("
                + std::to_string(weight.is_null()) + "/"
                + std::to_string(!found_order) +
                + ") from the XC energy density contraction "
                + expression->__str__()
            );
        }
        else {
            if (factors.empty()) {
                // We return nullptr when there is no generalized density
                // vectors, for example in the unperturbed case
                return std::make_tuple(
                    weight, order, SymEngine::RCP<const SymEngine::Basic>()
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
                    return std::make_tuple(weight, order, dens_vectors);
                }
                else {
                    throw SymEngine::SymEngineException(
                        "Invalid grid weights and/or XC functional derivatives "
                        + factors_left->__str__()
                        + " in "
                        + dens_vectors->__str__()
                    );
                }
            }
        }
    }

    std::pair<SymEngine::RCP<const OneElecOperator>, ExcContractionMap>
    extract_vxc_contraction(
        const SymEngine::RCP<const SymEngine::MatrixMul>& expression
    )
    {
        SymEngine::RCP<const OneElecOperator> Omega;
        ExcContractionMap energy_map;
        SymEngine::RCP<const SymEngine::Basic> scalar = SymEngine::one;
        for (const auto& arg: expression->get_args()) {
            // (Un)perturbed XC potential
            if (SymEngine::is_a_sub<const ExchCorrEnergy>(*arg)) {
                if (!energy_map.empty()) throw SymEngine::SymEngineException(
                    "Two XC potential contractions got from "
                    + expression->__str__()
                );
                auto op = SymEngine::rcp_dynamic_cast<const ExchCorrEnergy>(arg);
                energy_map = op->get_energy_map();
            }
            // Generalized overlap distribution
            else if (SymEngine::is_a_sub<const OneElecOperator>(*arg)) {
                if (!Omega.is_null()) throw SymEngine::SymEngineException(
                    "Two overlap distributions got from the XC potential contraction "
                    + expression->__str__()
                );
                Omega = SymEngine::rcp_dynamic_cast<const OneElecOperator>(arg);
            }
            else {
                scalar = SymEngine::mul(scalar, arg);
                //throw SymEngine::SymEngineException(
                //    "Invalid type from the XC potential contraction "
                //    + expression->__str__()
                //);
            }
        }
        // Consider the coefficient of matrix multiplication
        if (SymEngine::neq(*scalar, *SymEngine::one)) {
            for (auto& weight_map: energy_map) {
                for (auto& exc_map: weight_map.second)
                    exc_map.second = SymEngine::mul(exc_map.second, scalar);
            }
        }
        if (Omega.is_null() || energy_map.empty()) {
            throw SymEngine::SymEngineException(
                "Terms missing ("
                + std::to_string(Omega.is_null()) + "/"
                + std::to_string(energy_map.empty())
                + ") from the XC potential contraction "
                + expression->__str__()
            );
        }
        else {
            return std::make_pair(Omega, energy_map);
        }
    }

    void merge_exc_contraction(ExcContractionMap& map1, const ExcContractionMap& map2)
    {
        for (const auto& weight_map2: map2) {
            auto weight_map1 = map1.find(weight_map2.first);
            // New grid weight and corresponding contractions
            if (weight_map1 == map1.end()) {
                 map1.emplace(weight_map2);
            }
            else {
                for (const auto& exc_map2: weight_map2.second) {
                    auto exc_map1 = weight_map1->second.find(exc_map2.first);
                    // New XC energy derivative found
                    if (exc_map1 == weight_map1->second.end()) {
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

    bool eq_exc_contraction(
        const ExcContractionMap& map1, const ExcContractionMap& map2
    )
    {
        if (map1.size() == map2.size()) {
            auto weight_map1 = map1.begin();
            auto weight_map2 = map2.begin();
            for (; weight_map1 != map1.end(); ++weight_map1, ++weight_map2) {
                // Compare grid weights
                if (SymEngine::eq(*weight_map1->first, *weight_map2->first)) {
                    if (weight_map1->second.size() == weight_map2->second.size()) {
                        auto exc_map1 = weight_map1->second.begin();
                        auto exc_map2 = weight_map2->second.begin();
                        for (; exc_map1 != weight_map1->second.end();) {
                            if (exc_map1->first == exc_map2->first) {
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

    bool eq_vxc_contraction(
        const VxcContractionMap& map1, const VxcContractionMap& map2
    )
    {
        if (map1.size() == map2.size()) {
            auto energy_map1 = map1.begin();
            auto energy_map2 = map2.begin();
            for (; energy_map1 != map1.end(); ++energy_map1, ++energy_map2) {
                // Compare generalized overlap distribution
                if (SymEngine::eq(*energy_map1->first, *energy_map2->first)) {
                    if (!eq_exc_contraction(
                        energy_map1->second, energy_map2->second
                    )) return false;
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
