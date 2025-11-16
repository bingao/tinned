#include <stdio.h>
#include <string.h>

#include "tinned.h"
#include "tinned_cleanup.h"

int main(void) {
    TinnedErrorHandle_t* err = NULL;

    // Create perturbation a
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

    // Create perturbation b
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

    // Create dependencies with respect to perturbations
    PertMultichainHandle_t* dependencies = tinned_pert_multichain_new();
    bool ok = tinned_pert_multichain_insert(dependencies, pert_a, &err);
    if (!ok) {
        fprintf(
            stderr,
            "Failed to insert perturbation a into dependencies, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    ok = tinned_pert_multichain_insert(dependencies, pert_b, &err);
    if (!ok) {
        fprintf(
            stderr,
            "Failed to insert perturbation b into dependencies, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    // Create a one-electron operator
    ExprHandle_t* oper_1el = tinned_one_elec_operator_new("h", dependencies, &err);
    if (!oper_1el) {
        fprintf(
            stderr,
            "Failed to create one-electron operator, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    // Create a density matrix
    ExprHandle_t* ao_dens = tinned_wfn_parameter_new("D", &err);
    if (!ao_dens) {
        fprintf(
            stderr,
            "Failed to create density matrix, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    // Create two-electron energy
    ExprHandle_t* energy_2el = tinned_two_elec_energy_new("E_H", ao_dens, NULL, true, dependencies, NULL, &err);
    if (!energy_2el) {
        fprintf(
            stderr,
            "Failed to create two-electron energy, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    TINNED_SAFE_FREE_PERT_MULTICHAIN(dependencies);

    // Create the one-electron energy
    const ExprHandle_t* hD_terms[2] = {oper_1el, ao_dens};
    slice_ref_ExprHandle_const_ptr_t hD_slice = {
        .ptr = hD_terms,
        .len = 2,
    };
    ExprHandle_t* hD = tinned_matrix_mul_new(hD_slice, &err);
    if (!hD) {
        fprintf(
            stderr,
            "Failed to creat h*D, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    ExprHandle_t* energy_1el = tinned_trace_new(hD, &err);
    if (!energy_1el) {
        fprintf(
            stderr,
            "Failed to create one-electron energy, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    TINNED_SAFE_FREE_EXPR(oper_1el);
    TINNED_SAFE_FREE_EXPR(ao_dens);

    // Create the energy
    const ExprHandle_t* energy_terms[2] = {energy_1el, energy_2el};
    slice_ref_ExprHandle_const_ptr_t energy_slice = {
        .ptr = energy_terms,
        .len = 2,
    };
    ExprHandle_t* energy = tinned_add_new(energy_slice, &err);
    if (!energy) {
        fprintf(
            stderr,
            "Failed to create the energy, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    TINNED_SAFE_FREE_EXPR(energy_1el);
    TINNED_SAFE_FREE_EXPR(energy_2el);

    // Differentiate the energy and serialize it
    ExprHandle_t* energy_a = tinned_expr_differentiate(energy, pert_a, &err);
    if (!energy_a) {
        fprintf(
            stderr,
            "Failed to differentiate the energy w.r.t. perturbation a, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    char* str_energy = tinned_expr_serialize_json(energy, &err);
    if (!str_energy) {
        fprintf(
            stderr,
            "Failed to serialize the energy, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }
    char* str_energy_a = tinned_expr_serialize_json(energy_a, &err);
    if (!str_energy_a) {
        fprintf(
            stderr,
            "Failed to serialize the differentiated energy, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    fprintf(stdout, "E = %s\n", str_energy);
    fprintf(stdout, "E^{a} = %s\n", str_energy_a);

    // Cleanup
    TINNED_SAFE_FREE_STR(str_energy);
    TINNED_SAFE_FREE_STR(str_energy_a);
    TINNED_SAFE_FREE_EXPR(energy);
    TINNED_SAFE_FREE_EXPR(energy_a);
    TINNED_SAFE_FREE_PERTURBATION(pert_a);
    TINNED_SAFE_FREE_PERTURBATION(pert_b);
    if (err) TINNED_SAFE_FREE_ERR(err);

    return 0;
}
