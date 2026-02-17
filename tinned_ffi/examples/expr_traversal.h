#ifndef EXPR_TRAVERSAL_H
#define EXPR_TRAVERSAL_H

#include "tinned.h"

// Simple context structure passed through ctx
typedef struct {
    size_t depth;
    size_t node_count;
    size_t leaf_count;
} TraversalCtx;

//FIXME: Change name to assembling nodes
bool traversal_begin_node(void *ctx, ExprTag_t tag, size_t arity);
//FIXME: Leaf node should be some basic and non-assembling node
bool traversal_on_leaf(void *ctx, ExprTag_t tag, ExprHandle_t* expr);
bool traversal_end_node(void *ctx, ExprTag_t tag, size_t arity);

#endif