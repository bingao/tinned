#include <stdio.h>
#include <string.h>

#include "tinned.h"

int run_test_symbol(void) {
    TinnedErrorHandle_t* err = NULL;

    // Create Symbol("alpha")
    ExprHandle_t* sym = tinned_symbol_new("alpha", &err);
    if (!sym) {
        fprintf(stderr, "tinned_symbol_new failed\n");
        tinned_error_free(err);
        return 1;
    }

    // Get its name
    char* name = tinned_symbol_name(sym, &err);
    if (!name) {
        fprintf(stderr, "tinned_symbol_name failed\n");
        tinned_expr_free(sym);
        tinned_error_free(err);
        return 1;
    }

    // Verify
    int rc = 0;
    if (strcmp(name, "alpha") != 0) {
        fprintf(stderr, "unexpected symbol name: '%s'\n", name);
        rc = 1;
    }

    // Cleanup
    tinned_string_free(name);
    tinned_expr_free(sym);
    if (err) tinned_error_free(err);

    return rc;
}