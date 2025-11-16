typedef struct {
    FILE* fp;
    int count;
} PrintCtx;

static void on_item(void* ctx, uint32_t order, ExprHandle expr) {
    PrintCtx* C = (PrintCtx*)ctx;
    fprintf(C->fp, "order=%u  key=%s\n", order, tinned_expr_hash_key(&expr));
    C->count += 1;

    // You likely own the handle; free when done (per your API’s ownership rules).
    tinned_expr_free(&expr);
}

void example(ExprHandle* h, ExprHandle* s) {
    PrintCtx pc = { .fp = stdout, .count = 0 };
    SuperchainVisitor v = { .ctx = &pc, .on_item = on_item };
    TinnedErrorHandle_t* err = NULL;

    bool ok = tinned_expr_find_superchains_each(h, s, v, &err);
    if (!ok) { /* handle err */ }

    printf("Total items: %d\n", pc.count);
}

bool my_begin(void* ctx, ExprTag tag, size_t arity) {
    switch (tag) {
        case EXPR_TAG_ADD:
            /* ... */
            return true;
        case EXPR_TAG_ADJOINT_MAP:
            /* ... */
            return true;
        default:
            return false;
    }
}
