#include <stdio.h>
#include <string.h>

#include "tinned_cleanup.h"

int run_test_symbol(void) {
    TinnedErrorHandle_t* err = NULL;

    // Create Symbol("alpha")
    ExprHandle_t* symbol = tinned_symbol_new("alpha", &err);
    if (!symbol) {
        fprintf(stderr, "tinned_symbol_new failed\n");
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }

    // Get its name
    char* name = tinned_symbol_name(symbol, &err);
    if (!name) {
        fprintf(stderr, "tinned_symbol_name failed\n");
        TINNED_SAFE_FREE_EXPR(symbol);
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }

    // Verify
    if (strcmp(name, "alpha") != 0) {
        fprintf(stderr, "Unexpected symbol name: '%s'\n", name);
        TINNED_SAFE_FREE_EXPR(symbol);
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }
    TINNED_SAFE_FREE_STR(name);

    // Get its type name
    char* type_name = tinned_expr_type_name(symbol, &err);
    if (type_name) {
        fprintf(stdout, "Expr type is %s\n", type_name);
    } else {
        fprintf(stderr, "tinned_expr_type_name failed\n");
        TINNED_SAFE_FREE_EXPR(symbol);
        TINNED_SAFE_FREE_ERR(err);
        return 1;
    }
    TINNED_SAFE_FREE_STR(type_name);

    // Cleanup
    TINNED_SAFE_FREE_EXPR(symbol);
    if (err) TINNED_SAFE_FREE_ERR(err);

    return 0;
}
