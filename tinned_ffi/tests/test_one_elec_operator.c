#include <stdio.h>
#include <string.h>

#include "tinned_cleanup.h"

int run_test_one_elec_operator(void) {
    TinnedErrorHandle_t* err = NULL;

    // Create a perturbation
    const char* pert_name = "EL";
    ExprHandle_t* pert_freq = tinned_number_one_half();
    if (!pert_freq) {
        fprintf(stderr, "tinned_number_one_half failed\n");
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }
    PerturbationHandle_t* pert_el = tinned_perturbation_new(pert_name, pert_freq, &err);

    // Create dependencies with respect to perturbations
    PertMultichainHandle_t* dependencies = tinned_pert_multichain_new();
    bool ok = tinned_pert_multichain_insert(dependencies, pert_el, &err);
    if (!ok) {
        fprintf(stderr, "tinned_pert_multichain_insert failed\n");
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }
    TINNED_SAFE_FREE_EXPR(pert_freq);
    TINNED_SAFE_FREE_PERTURBATION(pert_el);

    // Create OneElecOperator("alpha")
    ExprHandle_t* op = tinned_one_elec_operator_new("alpha", dependencies, &err);
    if (!op) {
        fprintf(stderr, "tinned_one_elec_operator_new failed\n");
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }
    TINNED_SAFE_FREE_PERT_MULTICHAIN(dependencies);

    // Get the oeprator's name
    char* oper_name = tinned_one_elec_operator_name(op, &err);
    if (!oper_name) {
        fprintf(stderr, "tinned_one_elec_operator_name failed\n");
        TINNED_SAFE_FREE_EXPR(op);
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }

    // Verify
    if (strcmp(oper_name, "alpha") != 0) {
        fprintf(stderr, "Unexpected operator name: '%s'\n", oper_name);
        TINNED_SAFE_FREE_EXPR(op);
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }
    TINNED_SAFE_FREE_STR(oper_name);

    // Get the operator's type name
    char* oper_type = tinned_expr_type_name(op, &err);
    if (oper_type) {
        fprintf(stdout, "Expr type is %s\n", oper_type);
    } else {
        fprintf(stderr, "tinned_expr_type_name failed\n");
        TINNED_SAFE_FREE_EXPR(op);
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }
    TINNED_SAFE_FREE_STR(oper_type);

    // Cleanup
    TINNED_SAFE_FREE_EXPR(op);
    if (err) TINNED_SAFE_FREE_ERR(err);

    return 0;
}
