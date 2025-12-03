#include <stdio.h>
#include <string.h>

#include "tinned.h"
#include "tinned_cleanup.h"

#include "expr_traversal.h"

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
    PerturbationEntry_t pert_entries[2];
    pert_entries[0] = tinned_perturbation_entry_new(pert_a, 1);
    pert_entries[1] = tinned_perturbation_entry_new(pert_b, 2);
    PerturbationEntrySlice_t pert_slice = {
        .ptr = pert_entries,
        .len = 2,
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
    ExprHandle_t const * const hD_terms[2] = {oper_1el, ao_dens};
    ExprSlice_t hD_slice = {
        .ptr = hD_terms,
        .len = 2,
    };
    ExprHandle_t* hD = tinned_matrix_mul_new(&hD_slice, &err);
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
    ExprHandle_t const * const energy_terms[2] = {energy_1el, energy_2el};
    ExprSlice_t energy_slice = {
        .ptr = energy_terms,
        .len = 2,
    };
    ExprHandle_t* energy = tinned_add_new(&energy_slice, &err);
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

    fprintf(stdout, "E = %s\n\n", str_energy);
    fprintf(stdout, "E^{a} = %s\n\n", str_energy_a);

    // Traverse the expression tree
    TraversalCtx ctx = {
        .depth = 0,
        .node_count = 0,
        .leaf_count = 0,
    };

    CExprVisitor_t visitor = {
        .ctx       = &ctx,
        .begin_node = traversal_begin_node,
        .on_leaf    = traversal_on_leaf,
        .end_node   = traversal_end_node,
    };

    bool ok = tinned_walk_expr_postorder(energy_a, visitor, &err);
    if (!ok) {
        fprintf(
            stderr,
            "Failed to traverse the differentiated energy, with error message: %s\n",
            tinned_error_display(err)
        );
        return 1;
    }

    fprintf(stdout, "\nTraversal finished.\n");
    fprintf(stdout, "Visited nodes: %zu\n", ctx.node_count);
    fprintf(stdout, "Visited leaves: %zu\n", ctx.leaf_count);

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
