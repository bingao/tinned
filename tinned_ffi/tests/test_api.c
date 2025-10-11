#include <stdio.h>

#include "tinned.h"

typedef int (*test_fn)(void);
typedef struct { const char* name; test_fn fn; } test_entry;

/* prototypes from the list */
#define X(name) extern int name(void);
#include "all_tests.def"
#undef X

/* registry from the same list */
static const test_entry TESTS[] = {
#define X(name) { #name, name },
#include "all_tests.def"
#undef X
};

int main(void) {
    int fails = 0;
    size_t n = sizeof(TESTS) / sizeof(TESTS[0]);
    for (size_t i = 0; i < n; ++i) {
        int rc = TESTS[i].fn();
        if (rc != 0) {
            fprintf(stderr, "FAIL: %s (rc=%d)\n", TESTS[i].name, rc);
            fails++;
        }
    }
    if (fails == 0) { puts("All C FFI tests passed."); return 0; }
    fprintf(stderr, "%d test(s) failed.\n", fails);
    return 1;
}