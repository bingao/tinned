#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/PertDependency.hpp"
#include "Tinned/AdjointMap.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/ExpAdjointHamiltonian.hpp"

namespace Tinned
{
    ExpAdjointHamiltonian::ExpAdjointHamiltonian(
        const std::string& name,
        const SymEngine::RCP<const StateOperator>& state_operator,
        const SymEngine::RCP<const SymEngine::Basic>& hamiltonian
    ) : SymEngine::MatrixSymbol(name),
        state_operator_(state_operator)
    {
        if (SymEngine::is_a_sub<const AdjointMap>(*hamiltonian)) {
            auto adj_hamiltonain = SymEngine::rcp_dynamic_cast<const AdjointMap>(
                hamiltonian
            );
            if (adj_hamiltonain->size_x()>3) throw SymEngine::SymEngineException(
                "Exponential map does not allow fourfold or high-fold commutators"
            );
        }
        //FIXME: currently we only allow `OneElecOperator`
        else if (!SymEngine::is_a_sub<const OneElecOperator>(*hamiltonian)) {
            throw SymEngine::SymEngineException(
                "Exponential map currently only allows AdjointMap and OneElecOperator"
            );
        }
        hamiltonian_ = hamiltonian;
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ExpAdjointHamiltonian::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *state_operator_);
        SymEngine::hash_combine(seed, *hamiltonian_);
        return seed;
    }

    bool ExpAdjointHamiltonian::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const ExpAdjointHamiltonian>(o)) {
            auto& op = SymEngine::down_cast<const ExpAdjointHamiltonian&>(o);
            return get_name() == op.get_name()
                && state_operator_->__eq__(*op.state_operator_)
                && hamiltonian_->__eq__(*op.hamiltonian_);
        }
        return false;
    }

    int ExpAdjointHamiltonian::compare(const SymEngine::Basic &o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ExpAdjointHamiltonian>(o))
        auto& op = SymEngine::down_cast<const ExpAdjointHamiltonian&>(o);
        if (get_name() == op.get_name()) {
            int result = state_operator_->compare(*op.state_operator_);
            return result == 0 ? hamiltonian_->compare(*op.hamiltonian_) : result;
        }
        else {
            return get_name() < op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic ExpAdjointHamiltonian::get_args() const
    {
        return SymEngine::vec_basic({state_operator_, hamiltonian_});
    }

    SymEngine::RCP<const SymEngine::Basic> ExpAdjointHamiltonian::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        SymEngine::vec_basic terms;
        // Differentiate the Hamiltonian
        auto diff_hamiltonian = hamiltonian_->diff(s);
        if (!diff_hamiltonian->__eq__(*make_zero_operator())) {
            terms.push_back(SymEngine::make_rcp<const ExpAdjointHamiltonian>(
                get_name(), state_operator_, diff_hamiltonian
            ));
        }
        // Differentiate the state operator
        auto diff_oper = state_operator_->diff(s);
        if (!diff_oper->__eq__(*make_zero_operator())) {
            auto factors = SymEngine::vec_basic({SymEngine::minus_one});
            if (SymEngine::is_a_sub<const AdjointMap>(*hamiltonian_)) {
                auto adj_hamiltonain = SymEngine::rcp_dynamic_cast<const AdjointMap>(
                    hamiltonian_
                );
                if (adj_hamiltonain->size_x()<3) {
                    factors.push_back(SymEngine::make_rcp<const ExpAdjointHamiltonian>(
                        get_name(),
                        state_operator_,
                        SymEngine::make_rcp<const AdjointMap>(
                            *adj_hamiltonain, diff_oper
                        )
                    ));
                }
                // Fourfold commutators
                else {
                    factors.push_back(SymEngine::make_rcp<const AdjointMap>(
                        *adj_hamiltonain, diff_oper
                    ));
                }
            }
            // Here we do not check its type of `hamiltonian_` that has been
            // verified in the constructor
            else if (SymEngine::is_a_sub<const SymEngine::MatrixSymbol>(*hamiltonian_)) {
                auto name_hamiltonian = (
                    SymEngine::rcp_dynamic_cast<const SymEngine::MatrixSymbol>(hamiltonian_)
                )->get_name();
                factors.push_back(SymEngine::make_rcp<const ExpAdjointHamiltonian>(
                    get_name(),
                    state_operator_,
                    SymEngine::make_rcp<const AdjointMap>(
                        std::string("ad_{")
                            + state_operator_->get_name()
                            + std::string("}(")
                            + name_hamiltonian
                            + std::string(")"),
                        SymEngine::vec_basic({diff_oper}),
                        hamiltonian_
                    )
                ));
            }
            else {
                throw SymEngine::SymEngineException(
                    "ExpAdjointHamiltonian::diff_impl got an invalid Hamiltonian"
                );
            }
            terms.push_back(SymEngine::matrix_mul(factors));
        }
        if (terms.empty()) {
            return make_zero_operator();
        }
        else {
            return terms.size()==1 ? terms[0] : SymEngine::add(terms);
        }
    }
}
