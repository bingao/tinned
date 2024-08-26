#include "Tinned/TwoLevelAtom.hpp"

namespace Tinned
{
    TwoLevelOperator::TwoLevelOperator(
        const std::pair<SymEngine::RCP<const OneElecOperator>,
                        SymEngine::RCP<const SymEngine::MatrixExpr>>& H0,
        const std::map<SymEngine::RCP<const OneElecOperator>,
                       SymEngine::RCP<const SymEngine::MatrixExpr>,
                       SymEngine::RCPBasicKeyLess>& V,
        const std::pair<SymEngine::RCP<const OneElecOperator>,
                        SymEngine::RCP<const SymEngine::MatrixExpr>>& rho0
    ) : max_order_(7), H0_(H0), V_(V)
    {
        // Density matrix should be idempotent and has purity one
        auto Z = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
            SymEngine::matrix_add({
                SymEngine::matrix_mul({rho0.second, rho0.second}),
                SymEngine::matrix_mul({SymEngine::minus_one, rho0.second})
            })
        );
        for (const auto& val: Z->get_values())
            if (!is_zero_quantity(val)) throw SymEngine::SymEngineException(
                "Density matrix isn't idempotent: " + stringify(rho0.second)
            );
        if (!is_zero_quantity(SymEngine::add(SymEngine::minus_one, SymEngine::trace(rho0.second))))
            throw SymEngine::SymEngineException(
                "Density matrix doesn't have purity one: " + stringify(rho0.second)
            );
        rho0_ = rho0;
        SymEngine::vec_basic values;
        auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(H0.second);
        for (std::size_t i=0; i<op->nrows(); ++i)
            for (std::size_t j=0; j<op->ncols(); ++j) {
                if (i==j) {
                    values.push_back(SymEngine::zero);
                }
                else {
                    values.push_back(SymEngine::sub(op->get(i, i), op->get(j, j)));
                }
            }
        omega_ = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
            SymEngine::immutable_dense_matrix(op->nrows(), op->ncols(), values)
        );
        values.clear();
        for (const auto& oper: V_) {
            auto dependencies = oper.first->get_dependencies();
            if (dependencies.size()!=1) throw SymEngine::SymEngineException(
                "Each field operator should only depend on one perturbation: " + stringify(oper.first)
            );
            if (perturbations_.find(dependencies.begin()->first)!=perturbations_.end())
                throw SymEngine::SymEngineException(
                    "Field operators should depend on different perturbations: " + stringify(oper.first)
                );
            perturbations_.insert(dependencies.begin()->first);
            values.push_back(oper.second);
        }
        // All V's should commute
        for (std::size_t i=1; i<values.size(); ++i) {
            for (std::size_t j=0; j<i; ++j) {
                auto commutator = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
                    SymEngine::matrix_add({
                        SymEngine::matrix_mul({values[i], values[j]}),
                        SymEngine::matrix_mul({SymEngine::minus_one, values[j], values[i]})
                    })
                );
                for (const auto& val: commutator->get_values())
                    if (!is_zero_quantity(val)) throw SymEngine::SymEngineException(
                        "Field operators: " + stringify(values[i]) + " and "
                        + stringify(values[j]) + " do not commute"
                    );
            }
        }
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_hermitian_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A)
    {
        auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
        SymEngine::vec_basic values;
        for (std::size_t j=0; j<op->ncols(); ++j)
            for (std::size_t i=0; i<op->nrows(); ++i)
                values.push_back(
                    SymEngine::rcp_dynamic_cast<const SymEngine::Number>(op->get(i, j))->conjugate()
                );
        return SymEngine::immutable_dense_matrix(op->ncols(), op->nrows(), values);
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_1el_density(const OneElecDensity& x)
    {
        auto derivatives = x.get_derivatives();
        if (derivatives.size()>max_order_) SymEngine::SymEngineException(
            "Invalid order: " + std::to_string(derivatives.size())
            + "(" + std::to_string(max_order_) + ")"
        );
        switch (derivatives.size()) {
            case 0:
                return rho0_.second;
            default:
                if (derivatives.size()>rho_all_derivatives_.size()) {
                    for (unsigned int order=rho_all_derivatives_.size();
                         order<derivatives.size();
                         ++order) {
                        DensityDerivative rho_derivatives;
                        auto pert_permutation = PertPermutation(order+1, perturbations_);
                        bool remaining = true;
                        do {
                            auto permut_derivatives
                                = pert_permutation.get_derivatives(remaining);
                            auto freq_sum = get_frequency_sum(permut_derivatives);
                            // The last perturbation is the one for the
                            // external field, and we need to find the value of
                            // its corresponding operator
                            SymEngine::RCP<const SymEngine::MatrixExpr> oper_field;
                            bool not_found = true;
                            for (const auto& oper: V_) {
                                if (SymEngine::eq(
                                    *permut_derivatives.back(),
                                    *oper.first->get_dependencies().begin()->first
                                )) {
                                    oper_field = oper.second;
                                    not_found = false;
                                    break;
                                }
                            }
                            if (not_found) throw SymEngine::SymEngineException(
                                "Invalid perturbation for the external field "
                                + stringify(permut_derivatives.back())
                            );
                            // Check if the higher order derivatives of density
                            // matrix already computed partially
                            auto higher_derivatives = SymEngine::multiset_basic(
                                permut_derivatives.begin(), permut_derivatives.end()
                            );
                            int idx_higher = -1;
                            for (std::size_t i=0; i<rho_derivatives.size(); ++i)
                                if (SymEngine::unified_eq(rho_derivatives[i].first, higher_derivatives)) {
                                    idx_higher = i;
                                    break;
                                }
                            // Try to find the matched lower order derivatives
                            // of density matrix
                            not_found = true;
                            auto lower_derivatives = SymEngine::multiset_basic(
                                permut_derivatives.begin(), permut_derivatives.end()-1
                            );
                            for (const auto& rho_lower: rho_all_derivatives_[order]) {
                                if (SymEngine::unified_eq(rho_lower.first, lower_derivatives)) {
                                    auto rho_commutator
                                        = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(
                                              SymEngine::matrix_add({
                                                  SymEngine::matrix_mul({
                                                      oper_field, rho_lower.second
                                                  }),
                                                  SymEngine::matrix_mul({
                                                      SymEngine::minus_one,
                                                      rho_lower.second,
                                                      oper_field
                                                  })
                                              })
                                          );
                                    SymEngine::vec_basic val_rho;
                                    for (std::size_t i=0; i<omega_->nrows(); ++i)
                                        for (std::size_t j=0; j<omega_->ncols(); ++j)
                                            val_rho.push_back(
                                                SymEngine::div(
                                                    rho_commutator->get(i, j),
                                                    SymEngine::sub(
                                                        freq_sum,
                                                        omega_->get(i, j)
                                                    )
                                                )
                                            );
                                    if (idx_higher>=0) {
                                        rho_derivatives[idx_higher].second
                                            = SymEngine::matrix_add({
                                                  rho_derivatives[idx_higher].second,
                                                  SymEngine::immutable_dense_matrix(
                                                      omega_->nrows(),
                                                      omega_->ncols(),
                                                      val_rho
                                                  )
                                              });
                                    }
                                    else {
                                        rho_derivatives.push_back(
                                            std::make_pair(
                                                higher_derivatives,
                                                SymEngine::immutable_dense_matrix(
                                                    omega_->nrows(),
                                                    omega_->ncols(),
                                                    val_rho
                                                )
                                            )
                                        );
                                    }
                                    not_found = false;
                                    break;
                                }
                            }
                            if (not_found) {
                                std::string str_lower_derivatives;
                                for (const auto& p: permut_derivatives)
                                   str_lower_derivatives += stringify(p) + ",";
                                throw SymEngine::SymEngineException(
                                    "Invalid lower order derivatives: "
                                    + str_lower_derivatives
                                    + " of density matrix: "
                                    + stringify(x)
                                );
                            }
                        } while (remaining);
                        rho_all_derivatives_[order+1] = rho_derivatives;
                    }
                }
                for (const auto& rho_derivatives: rho_all_derivatives_[derivatives.size()])
                    if (SymEngine::unified_eq(derivatives, rho_derivatives.first))
                        return rho_derivatives.second;
                throw SymEngine::SymEngineException(
                    "Density matrix with invalid derivatives: " + stringify(x)
                );
        }
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_1el_operator(const OneElecOperator& x)
    {
        if (x.get_name()==H0_.first->get_name()) {
            auto derivatives = x.get_derivatives();
            if (derivatives.empty()) {
                return H0_.second;
            }
            else {
                return make_zero_matrix();
            }
        }
        else {
            for (const auto& oper: V_) {
                if (x.get_name()==oper.first->get_name()) {
                    auto derivatives = x.get_derivatives();
                    if (derivatives.empty()) {
                        return SymEngine::matrix_mul({
                            oper.first->get_dependencies().begin()->first,
                            oper.second
                        });
                    }
                    else if (derivatives.size()==1) {
                        return oper.second;
                    }
                    else {
                        return make_zero_matrix();
                    }
                }
            }
        }
        throw SymEngine::SymEngineException("Invalid operator " + stringify(x));
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_temporum_operator(const TemporumOperator& x)
    {
        result_ = apply(x.get_target());
        eval_oper_scale(x.get_frequency(), result_);
        return result_;
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_conjugate_matrix(const SymEngine::RCP<const SymEngine::MatrixExpr>& A)
    {
        auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
        SymEngine::vec_basic values;
        for (std::size_t i=0; i<op->nrows(); ++i)
            for (std::size_t j=0; j<op->ncols(); ++j)
                values.push_back(
                    SymEngine::rcp_dynamic_cast<const SymEngine::Number>(op->get(i, j))->conjugate()
                );
        return SymEngine::immutable_dense_matrix(op->nrows(), op->ncols(), values);
    }

    SymEngine::RCP<const SymEngine::MatrixExpr>
    TwoLevelOperator::eval_transpose(const SymEngine::RCP<const SymEngine::MatrixExpr>& A)
    {
        auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
        SymEngine::vec_basic values;
        for (std::size_t j=0; j<op->ncols(); ++j)
            for (std::size_t i=0; i<op->nrows(); ++i)
                values.push_back(op->get(i, j));
        return SymEngine::immutable_dense_matrix(op->ncols(), op->nrows(), values);
    }

    void TwoLevelOperator::eval_oper_addition(
        SymEngine::RCP<const SymEngine::MatrixExpr>& A,
        const SymEngine::RCP<const SymEngine::MatrixExpr>& B
    )
    {
        auto op_A = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
        auto op_B = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(B);
        SYMENGINE_ASSERT(
            op_A->nrows()==op_B->nrows() && op_A->ncols()==op_B->ncols()
        )
        SymEngine::vec_basic values;
        for (std::size_t i=0; i<op_A->nrows(); ++i)
            for (std::size_t j=0; j<op_A->ncols(); ++j)
                values.push_back(
                    SymEngine::add(op_A->get(i, j), op_B->get(i, j))
                );
        A = SymEngine::immutable_dense_matrix(op_A->nrows(), op_A->ncols(), values);
    }

    SymEngine::RCP<const SymEngine::MatrixExpr> TwoLevelOperator::eval_oper_multiplication(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& A,
        const SymEngine::RCP<const SymEngine::MatrixExpr>& B
    )
    {
        return SymEngine::matrix_mul({A, B});
    }

    void TwoLevelOperator::eval_oper_scale(
        const SymEngine::RCP<const SymEngine::Number>& scalar,
        SymEngine::RCP<const SymEngine::MatrixExpr>& A
    )
    {
        auto op = SymEngine::rcp_dynamic_cast<const SymEngine::ImmutableDenseMatrix>(A);
        SymEngine::vec_basic values;
        for (std::size_t i=0; i<op->nrows(); ++i)
            for (std::size_t j=0; j<op->ncols(); ++j)
                values.push_back(SymEngine::mul(scalar, op->get(i, j)));
        A = SymEngine::immutable_dense_matrix(op->nrows(), op->ncols(), values);
    }

    TwoLevelFunction::TwoLevelFunction(
        const std::pair<SymEngine::RCP<const OneElecOperator>,
                        SymEngine::RCP<const SymEngine::MatrixExpr>>& H0,
        const std::map<SymEngine::RCP<const OneElecOperator>,
                       SymEngine::RCP<const SymEngine::MatrixExpr>,
                       SymEngine::RCPBasicKeyLess>& V,
        const std::pair<SymEngine::RCP<const OneElecOperator>,
                        SymEngine::RCP<const SymEngine::MatrixExpr>>& rho0
    ) : FunctionEvaluator<SymEngine::RCP<const SymEngine::Basic>,
                          SymEngine::RCP<const SymEngine::MatrixExpr>>(
            std::make_shared<TwoLevelOperator>(TwoLevelOperator(H0, V, rho0))
        ) {}

    SymEngine::RCP<const SymEngine::Basic> TwoLevelFunction::eval_trace(
        const SymEngine::RCP<const SymEngine::MatrixExpr>& A
    )
    {
        return SymEngine::trace(A);
    }

    void TwoLevelFunction::eval_fun_addition(
        SymEngine::RCP<const SymEngine::Basic>& f,
        const SymEngine::RCP<const SymEngine::Basic>& g
    )
    {
        f = SymEngine::add(f, g);
    }

    void TwoLevelFunction::eval_fun_scale(
        const SymEngine::RCP<const SymEngine::Number>& scalar,
        SymEngine::RCP<const SymEngine::Basic>& f
    )
    {
        f = SymEngine::mul(scalar, f);
    }
}
