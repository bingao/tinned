#include <stdio.h>
#include <string.h>

#include "tinned.h"
#include "tinned_cleanup.h"

void eval_xc_energy(const ExprHandle_t* energy)
{
    TinnedErrorHandle_t* err = NULL;

    ExprHandle_t* xc_energy = tinned_exch_corr_energy_xc_energy(energy, &err);
    if (!xc_energy) {
        fprintf(
            stderr,
            "Failed to get exchange-correlation energy functional, with error message: %s\n",
            tinned_error_display(err)
        );
        return;
    }

    char* xc_energy_type = tinned_expr_type_name(xc_energy, &err);
    if (!xc_energy_type) {
        fprintf(
            stderr,
            "Failed to get type of the exchange-correlation energy functional, with error message: %s\n",
            tinned_error_display(err)
        );
        return;
    }

    if (strcmp(xc_energy_type, "Mul") == 0) {
        char* str_xc_energy = tinned_expr_serialize_json(xc_energy, &err);
        if (!str_xc_energy) {
            fprintf(
                stderr,
                "Failed to serialize the exchange-correlation energy functional, with error message: %s\n",
                tinned_error_display(err)
            );
            return;
        }
        fprintf(stdout, "Exc = %s\n\n", str_xc_energy);
        TINNED_SAFE_FREE_STR(str_xc_energy);
    } else {
        Vec_ExprHandle_ptr_t xc_terms = tinned_add_terms(xc_energy, &err);
        if (err != NULL) {
            fprintf(
                stderr,
                "Failed to get terms of the exchange-correlation energy functional, with error message: %s\n",
                tinned_error_display(err)
            );
            return;
        }
        // Iterate over the terms
        fprintf(stdout, "Exc has %ld terms\n\n", xc_terms.len);
        for (size_t i = 0; i < xc_terms.len; ++i) {
            ExprHandle_t *xc_term = xc_terms.ptr[i];
            char* str_xc_term = tinned_expr_serialize_json(xc_term, &err);
            if (!str_xc_term) {
                fprintf(
                    stderr,
                    "Failed to serialize a term of the exchange-correlation energy functional, with error message: %s\n",
                    tinned_error_display(err)
                );
                return;
            }
            fprintf(stdout, "Exc[%ld] = %s\n\n", i, str_xc_term);
            TINNED_SAFE_FREE_STR(str_xc_term);
        }
        tinned_expr_vec_free(xc_terms);
    }

    TINNED_SAFE_FREE_STR(xc_energy_type);
    if (err) TINNED_SAFE_FREE_ERR(err);
}

void eval_xc_potential(const ExprHandle_t* potential)
{
    TinnedErrorHandle_t* err = NULL;

    ExprHandle_t* xc_potential = tinned_exch_corr_potential_xc_potential(potential, &err);
    if (!xc_potential) {
        fprintf(
            stderr,
            "Failed to get exchange-correlation potential operator, with error message: %s\n",
            tinned_error_display(err)
        );
        return;
    }

    char* xc_potential_type = tinned_expr_type_name(xc_potential, &err);
    if (!xc_potential_type) {
        fprintf(
            stderr,
            "Failed to get type of the exchange-correlation potential operator, with error message: %s\n",
            tinned_error_display(err)
        );
        return;
    }

    if (strcmp(xc_potential_type, "MatrixMul") == 0) {
        char* str_xc_potential = tinned_expr_serialize_json(xc_potential, &err);
        if (!str_xc_potential) {
            fprintf(
                stderr,
                "Failed to serialize the exchange-correlation potential operator, with error message: %s\n",
                tinned_error_display(err)
            );
            return;
        }
        fprintf(stdout, "Vxc = %s\n\n", str_xc_potential);
        TINNED_SAFE_FREE_STR(str_xc_potential);
    } else {
        Vec_ExprHandle_ptr_t xc_terms = tinned_matrix_add_terms(xc_potential, &err);
        if (err != NULL) {
            fprintf(
                stderr,
                "Failed to get terms of the exchange-correlation potential operator, with error message: %s\n",
                tinned_error_display(err)
            );
            return;
        }
        // Iterate over the terms
        fprintf(stdout, "Vxc has %ld terms\n\n", xc_terms.len);
        for (size_t i = 0; i < xc_terms.len; ++i) {
            ExprHandle_t *xc_term = xc_terms.ptr[i];
            char* str_xc_term = tinned_expr_serialize_json(xc_term, &err);
            if (!str_xc_term) {
                fprintf(
                    stderr,
                    "Failed to serialize a term of the exchange-correlation potential operator, with error message: %s\n",
                    tinned_error_display(err)
                );
                return;
            }
            fprintf(stdout, "Vxc[%ld] = %s\n\n", i, str_xc_term);
            TINNED_SAFE_FREE_STR(str_xc_term);
        }
        tinned_expr_vec_free(xc_terms);
    }

    TINNED_SAFE_FREE_STR(xc_potential_type);
    if (err) TINNED_SAFE_FREE_ERR(err);
}

// This is just a simple illustration, probably not for pratical use
int main(void)
{
    TinnedErrorHandle_t* err = NULL;

    // Make perturbations
    ExprHandle_t* freq_a = tinned_symbol_new("omega_a", &err);
    if (!freq_a) {
        fprintf(
            stderr,
            "Failed to create frequency a, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    PerturbationHandle_t* pert_a = tinned_perturbation_new("a", freq_a, &err);
    if (!pert_a) {
        fprintf(
            stderr,
            "Failed to create perturbation a, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    TINNED_SAFE_FREE_EXPR(freq_a);

    ExprHandle_t* freq_b = tinned_symbol_new("omega_b", &err);
    if (!freq_b) {
        fprintf(
            stderr,
            "Failed to create frequency b, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    PerturbationHandle_t* pert_b = tinned_perturbation_new("b", freq_b, &err);
    if (!pert_b) {
        fprintf(
            stderr,
            "Failed to create perturbation b, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    TINNED_SAFE_FREE_EXPR(freq_b);

    ExprHandle_t* freq_c = tinned_symbol_new("omega_c", &err);
    if (!freq_c) {
        fprintf(
            stderr,
            "Failed to create frequency c, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    PerturbationHandle_t* pert_c = tinned_perturbation_new("c", freq_c, &err);
    if (!pert_c) {
        fprintf(
            stderr,
            "Failed to create perturbation c, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    TINNED_SAFE_FREE_EXPR(freq_c);

    ExprHandle_t* freq_d = tinned_symbol_new("omega_d", &err);
    if (!freq_d) {
        fprintf(
            stderr,
            "Failed to create frequency d, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    PerturbationHandle_t* pert_d = tinned_perturbation_new("d", freq_d, &err);
    if (!pert_d) {
        fprintf(
            stderr,
            "Failed to create perturbation d, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    TINNED_SAFE_FREE_EXPR(freq_d);

    // Make perturbations' dependecy and their maximum differentiated orders.
    // You can make different dependecies for different operators.
    PerturbationEntry_t pert_entries[4];
    pert_entries[0] = tinned_perturbation_entry_new(pert_a, 99);
    pert_entries[1] = tinned_perturbation_entry_new(pert_b, 99);
    pert_entries[2] = tinned_perturbation_entry_new(pert_c, 99);
    pert_entries[3] = tinned_perturbation_entry_new(pert_d, 99);
    PerturbationEntrySlice_t pert_slice = {
        .ptr = pert_entries,
        .len = 4,
    };

    PertMultichainHandle_t* dependencies = tinned_pert_multichain_from_entries(pert_slice, &err);
    if (!dependencies) {
        fprintf(
            stderr,
            "Failed to create dependencies, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    // Make one-electron spin-orbital density matrix
    ExprHandle_t* D = tinned_wfn_parameter_new("D", &err);
    if (!D) {
        fprintf(
            stderr,
            "Failed to create density matrix, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    // Make grid weight
    ExprHandle_t* weight = tinned_non_elec_function_new("weight", dependencies, &err);
    if (!weight) {
        fprintf(
            stderr,
            "Failed to create grid weight, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    // Make generalized overlap distribution
    ExprHandle_t* Omega = tinned_one_elec_operator_new("Omega", dependencies, &err);
    if (!Omega) {
        fprintf(
            stderr,
            "Failed to create generalized overlap distribution, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    // Make exchange-correlation energy functional and potential operator
    ExprHandle_t* Exc = tinned_exch_corr_energy_new("GGA", weight, D, Omega, &err);
    if (!Exc) {
        fprintf(
            stderr,
            "Failed to create exchange-correlation energy functional, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    ExprHandle_t* Vxc = tinned_exch_corr_potential_new("GGA", weight, D, Omega, &err);
    if (!Vxc) {
        fprintf(
            stderr,
            "Failed to create exchange-correlation potential operator, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    TINNED_SAFE_FREE_EXPR(D);
    TINNED_SAFE_FREE_EXPR(weight);
    TINNED_SAFE_FREE_EXPR(Omega);

    // (0) Unperturbed case
    eval_xc_energy(Exc);
    eval_xc_potential(Vxc);

    // (1) The first order XC energy density derivative
    ExprHandle_t* Exc_a = tinned_expr_differentiate(Exc, pert_a, &err);
    if (!Exc_a) {
        fprintf(
            stderr,
            "Failed to compute Exc^a, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    ExprHandle_t* Vxc_a = tinned_expr_differentiate(Vxc, pert_a, &err);
    if (!Vxc_a) {
        fprintf(
            stderr,
            "Failed to compute Vxc^a, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    eval_xc_energy(Exc_a);
    eval_xc_potential(Vxc_a);

    TINNED_SAFE_FREE_EXPR(Exc);
    TINNED_SAFE_FREE_EXPR(Vxc);

    // (2) The second order XC energy density derivative
    ExprHandle_t* Exc_ab = tinned_expr_differentiate(Exc_a, pert_b, &err);
    if (!Exc_ab) {
        fprintf(
            stderr,
            "Failed to compute Exc^ab, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    ExprHandle_t* Vxc_ab = tinned_expr_differentiate(Vxc_a, pert_b, &err);
    if (!Vxc_ab) {
        fprintf(
            stderr,
            "Failed to compute Vxc^ab, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    eval_xc_energy(Exc_ab);
    eval_xc_potential(Vxc_ab);

    TINNED_SAFE_FREE_EXPR(Exc_a);
    TINNED_SAFE_FREE_EXPR(Vxc_a);

    // (3) The third order XC energy density derivative
    ExprHandle_t* Exc_abc = tinned_expr_differentiate(Exc_ab, pert_c, &err);
    if (!Exc_abc) {
        fprintf(
            stderr,
            "Failed to compute Exc^abc, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    ExprHandle_t* Vxc_abc = tinned_expr_differentiate(Vxc_ab, pert_c, &err);
    if (!Vxc_abc) {
        fprintf(
            stderr,
            "Failed to compute Vxc^abc, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    eval_xc_energy(Exc_abc);
    eval_xc_potential(Vxc_abc);

    TINNED_SAFE_FREE_EXPR(Exc_ab);
    TINNED_SAFE_FREE_EXPR(Vxc_ab);

    // (4) The fourth order XC energy density derivative
    ExprHandle_t* Exc_abcd = tinned_expr_differentiate(Exc_abc, pert_d, &err);
    if (!Exc_abcd) {
        fprintf(
            stderr,
            "Failed to compute Exc^abcd, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    ExprHandle_t* Vxc_abcd = tinned_expr_differentiate(Vxc_abc, pert_d, &err);
    if (!Vxc_abcd) {
        fprintf(
            stderr,
            "Failed to compute Vxc^abcd, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    eval_xc_energy(Exc_abcd);
    eval_xc_potential(Vxc_abcd);

    TINNED_SAFE_FREE_EXPR(Exc_abc);
    TINNED_SAFE_FREE_EXPR(Vxc_abc);
    TINNED_SAFE_FREE_EXPR(Exc_abcd);
    TINNED_SAFE_FREE_EXPR(Vxc_abcd);

    TINNED_SAFE_FREE_PERTURBATION(pert_a);
    TINNED_SAFE_FREE_PERTURBATION(pert_b);
    TINNED_SAFE_FREE_PERTURBATION(pert_c);
    TINNED_SAFE_FREE_PERTURBATION(pert_d);
    if (err) TINNED_SAFE_FREE_ERR(err);

    return 0;
}
