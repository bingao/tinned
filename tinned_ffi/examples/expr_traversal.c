#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "tinned_cleanup.h"
#include "expr_traversal.h"

// Optional helper for indentation
static void print_indent(size_t depth) {
    for (size_t i = 0; i < depth; ++i) {
        fputs("    ", stdout);
    }
}

// ctx: pointer to TraversalCtx
// tag: expression tag for the current node
// arity: number of children (0 for leaf, >0 for non-leaf)
bool traversal_begin_node(void *ctx, ExprTag_t tag, size_t arity) {
    TraversalCtx *state = (TraversalCtx *)ctx;

    print_indent(state->depth);

    // Only (matrix) addition and multiplication are allowed here
    switch (tag)
    {
        case EXPR_TAG_ADD:
            printf("traversal_begin_node() gets an addition with %zu terms\n", arity);
            break;
        case EXPR_TAG_MATRIX_ADD:
            printf("traversal_begin_node() gets a matrix addition with %zu terms\n", arity);
            break;
        case EXPR_TAG_MATRIX_MUL:
            printf("traversal_begin_node() gets a matrix multiplication with %zu factors\n", arity);
            break;
        case EXPR_TAG_MUL:
            printf("traversal_begin_node() gets a multiplication with %zu factors\n", arity);
            break;
        default:
            return false;
    }

    state->node_count += 1;
    state->depth += 1;

    return true;
}

// ctx: pointer to TraversalCtx
// tag: expression tag for the leaf
// expr: owned ExprBox, must be freed here with `TINNED_SAFE_FREE_EXPR`
bool traversal_on_leaf(void *ctx, ExprTag_t tag, ExprHandle_t* expr) {
    TraversalCtx *state = (TraversalCtx *)ctx;

    print_indent(state->depth);

    switch (tag)
    {
        case EXPR_TAG_ADJOINT_MAP:
            // Here, one can access all functions for adjoint map from Tinned
            printf("traversal_on_leaf() gets an adjoint map\n");
            break;
        case EXPR_TAG_COMPOSITION:
            printf("traversal_on_leaf() gets a composition\n");
            break;
        case EXPR_TAG_CONJUGATE:
            printf("traversal_on_leaf() gets a conjugate\n");
            break;
        case EXPR_TAG_DOT_PRODUCT:
            printf("traversal_on_leaf() gets a dot product\n");
            break;
        case EXPR_TAG_EXCH_CORR_ENERGY:
            printf("traversal_on_leaf() gets an exchange-correlation energy\n");
            break;
        case EXPR_TAG_EXCH_CORR_POTENTIAL:
            printf("traversal_on_leaf() gets an exchange-correlation potential\n");
            break;
        case EXPR_TAG_EXP_ADJOINT_MAP:
            printf("traversal_on_leaf() gets an exponential adjoint map\n");
            break;
        case EXPR_TAG_HERMITIAN_TRANSPOSE:
            printf("traversal_on_leaf() gets a Hermitian transpose\n");
            break;
        case EXPR_TAG_LAG_MULTIPLIER:
            printf("traversal_on_leaf() gets a Lagrangian multiplier\n");
            break;
        case EXPR_TAG_NON_ELEC_FUNCTION:
            printf("traversal_on_leaf() gets a non-electron function\n");
            break;
        case EXPR_TAG_NUMBER: {
            TinnedErrorHandle_t* err = NULL;
            char* str_number = tinned_expr_display(expr, &err);
            if (!str_number) {
                printf(
                    "Failed to serialize a number with error message %s\n",
                    tinned_error_display(err)
                );
                return false;
            }
            printf("traversal_on_leaf() gets a number %s\n", str_number);
            TINNED_SAFE_FREE_STR(str_number);
            break;
        }
        case EXPR_TAG_ONE_ELEC_MATRIX:
            printf("traversal_on_leaf() gets a one-electron matrix\n");
            break;
        case EXPR_TAG_POWER:
            printf("traversal_on_leaf() gets a power\n");
            break;
        case EXPR_TAG_RESIDUE_PARAMETER:
            printf("traversal_on_leaf() gets a residue parameter\n");
            break;
        case EXPR_TAG_SYMBOL: {
            TinnedErrorHandle_t* err = NULL;
            char* name = tinned_symbol_name(expr, &err);
            if (!name) {
                printf(
                    "Failed to get the name of a symbol with error message %s\n",
                    tinned_error_display(err)
                );
                return false;
            }
            printf("traversal_on_leaf() gets a symbol %s\n", name);
            TINNED_SAFE_FREE_STR(name);
            break;
        }
        case EXPR_TAG_TIME_EVOLUTION:
            printf("traversal_on_leaf() gets a time-evolution operator\n");
            break;
        case EXPR_TAG_BASIS_TIME_EVOLUTION:
            printf("traversal_on_leaf() gets a time-evolution overlap\n");
            break;
        case EXPR_TAG_TRACE: {
            TinnedErrorHandle_t* err = NULL;
            ExprHandle_t* argument = tinned_trace_argument(expr, &err);
            if (!argument) {
                printf(
                    "Failed to get the argument of a trace with error message %s\n",
                    tinned_error_display(err)
                );
                return false;
            }
            char* str_argument = tinned_expr_display(argument, &err);
            if (!str_argument) {
                printf(
                    "Failed to serialize the argument with error message %s\n",
                    tinned_error_display(err)
                );
                return false;
            }
            printf("traversal_on_leaf() gets a trace with the argument %s\n", str_argument);
            TINNED_SAFE_FREE_EXPR(argument);
            TINNED_SAFE_FREE_STR(str_argument);
            break;
        }
        case EXPR_TAG_TRANSPOSE:
            printf("traversal_on_leaf() gets a transpose\n");
            break;
        case EXPR_TAG_AO_TWO_ELEC_ENERGY:
            printf("traversal_on_leaf() gets a two-electron energy\n");
            break;
        case EXPR_TAG_AO_TWO_ELEC_MATRIX:
            printf("traversal_on_leaf() gets a two-electron matrix\n");
            break;
        case EXPR_TAG_WFN_PARAMETER:
            printf("traversal_on_leaf() gets a wave function parameter\n");
            break;
        case EXPR_TAG_ZERO_OPERATOR:
            printf("traversal_on_leaf() gets a zero operator\n");
            break;
        default:
            return false;
    }

    state->leaf_count += 1;

    // We own expr and must free it here.
    TINNED_SAFE_FREE_EXPR(expr);

    return true;
}

// ctx: pointer to TraversalCtx
// tag: expression tag
// arity: same arity as in traversal_begin_node()
bool traversal_end_node(void *ctx, ExprTag_t tag, size_t arity) {
    TraversalCtx *state = (TraversalCtx *)ctx;

    // Depth was already incremented in traversal_begin_node(), so decrement here
    if (state->depth > 0) {
        state->depth -= 1;
    }

    print_indent(state->depth);

    switch (tag)
    {
        case EXPR_TAG_ADD:
            printf("traversal_end_node() gets an addition with %zu terms\n", arity);
            break;
        case EXPR_TAG_MATRIX_ADD:
            printf("traversal_end_node() gets a matrix addition with %zu terms\n", arity);
            break;
        case EXPR_TAG_MATRIX_MUL:
            printf("traversal_end_node() gets a matrix multiplication with %zu factors\n", arity);
            break;
        case EXPR_TAG_MUL:
            printf("traversal_end_node() gets a multiplication with %zu factors\n", arity);
            break;
        default:
            return false;
    }

    return true;
}
