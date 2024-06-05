#include <string>

#include <symengine/add.h>
#include <symengine/constants.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/symengine_assert.h>
#include <symengine/symengine_casts.h>
#include <symengine/symengine_exception.h>

#include "Tinned/AdjointMap.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/ClusterConjHamiltonian.hpp"
#include "Tinned/ZeroOperator.hpp"

namespace Tinned
{
    ClusterConjHamiltonian::ClusterConjHamiltonian(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& cluster_operator,
        const SymEngine::RCP<const SymEngine::Basic>& hamiltonian
    ) : SymEngine::MatrixSymbol(std::string("e^{ad}")),
        cluster_operator_(cluster_operator)
    {
        if (SymEngine::is_a_sub<const AdjointMap>(*hamiltonian)) {
            auto ad_hamiltonain = SymEngine::rcp_dynamic_cast<const AdjointMap>(
                hamiltonian
            );
            if (ad_hamiltonain->size_x()>3) throw SymEngine::SymEngineException(
                "ClusterConjHamiltonian does not allow fourfold or high-fold commutators"
            );
        }
        //FIXME: currently we allow only `OneElecOperator`
        else if (!SymEngine::is_a_sub<const OneElecOperator>(*hamiltonian)) {
            throw SymEngine::SymEngineException(
                "ClusterConjHamiltonian currently allows only AdjointMap and OneElecOperator"
            );
        }
        hamiltonian_ = hamiltonian;
        SYMENGINE_ASSIGN_TYPEID()
    }

    SymEngine::hash_t ClusterConjHamiltonian::__hash__() const
    {
        SymEngine::hash_t seed = SymEngine::MatrixSymbol::__hash__();
        SymEngine::hash_combine(seed, *cluster_operator_);
        SymEngine::hash_combine(seed, *hamiltonian_);
        return seed;
    }

    bool ClusterConjHamiltonian::__eq__(const SymEngine::Basic& o) const
    {
        if (SymEngine::is_a_sub<const ClusterConjHamiltonian>(o)) {
            auto& op = SymEngine::down_cast<const ClusterConjHamiltonian&>(o);
            return get_name()==op.get_name()
                && cluster_operator_->__eq__(*op.cluster_operator_)
                && hamiltonian_->__eq__(*op.hamiltonian_);
        }
        return false;
    }

    int ClusterConjHamiltonian::compare(const SymEngine::Basic& o) const
    {
        SYMENGINE_ASSERT(SymEngine::is_a_sub<const ClusterConjHamiltonian>(o))
        auto& op = SymEngine::down_cast<const ClusterConjHamiltonian&>(o);
        if (get_name()==op.get_name()) {
            int result = cluster_operator_->compare(*op.cluster_operator_);
            return result==0 ? hamiltonian_->compare(*op.hamiltonian_) : result;
        }
        else {
            return get_name()<op.get_name() ? -1 : 1;
        }
    }

    SymEngine::vec_basic ClusterConjHamiltonian::get_args() const
    {
        return SymEngine::vec_basic({cluster_operator_, hamiltonian_});
    }

    SymEngine::RCP<const SymEngine::Basic> ClusterConjHamiltonian::diff_impl(
        const SymEngine::RCP<const SymEngine::Symbol>& s
    ) const
    {
        SymEngine::vec_basic terms;
        // Differentiate the Hamiltonian
        auto diff_hamiltonian = hamiltonian_->diff(s);
        if (!diff_hamiltonian->__eq__(*make_zero_operator())) {
            terms.push_back(SymEngine::make_rcp<const ClusterConjHamiltonian>(
                cluster_operator_, diff_hamiltonian
            ));
        }
        // Differentiate the cluster operator
        auto diff_oper = cluster_operator_->diff(s);
        if (!diff_oper->__eq__(*make_zero_operator())) {
            auto factors = SymEngine::vec_basic({SymEngine::minus_one});
            if (SymEngine::is_a_sub<const AdjointMap>(*hamiltonian_)) {
                auto ad_hamiltonain = SymEngine::rcp_dynamic_cast<const AdjointMap>(
                    hamiltonian_
                );
                if (ad_hamiltonain->size_x()<3) {
                    factors.push_back(SymEngine::make_rcp<const ClusterConjHamiltonian>(
                        cluster_operator_,
                        SymEngine::make_rcp<const AdjointMap>(*ad_hamiltonain, diff_oper)
                    ));
                }
                // Fourfold commutators
                else {
                    factors.push_back(SymEngine::make_rcp<const AdjointMap>(
                        *ad_hamiltonain, diff_oper
                    ));
                }
            }
            // Here we do not check the type of `hamiltonian_`, which has been
            // verified in the constructor
            else if (SymEngine::is_a_sub<const SymEngine::MatrixSymbol>(*hamiltonian_)) {
                factors.push_back(SymEngine::make_rcp<const ClusterConjHamiltonian>(
                    cluster_operator_,
                    SymEngine::make_rcp<const AdjointMap>(
                        SymEngine::vec_basic({diff_oper}), hamiltonian_
                    )
                ));
            }
            else {
                throw SymEngine::SymEngineException(
                    "ClusterConjHamiltonian::diff_impl got an invalid Hamiltonian"
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
