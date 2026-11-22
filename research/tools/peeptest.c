/* Upper-bound probe: apply exact windowed peephole (all 11 ops) with window K to an existing op file.
 * usage: peeptest K OPSFILE NUMBERS...   prints "before after" */
#include "../c/ps.h"
int main(int argc, char **argv) {
    int K = atoi(argv[1]); FILE *f = fopen(argv[2], "r"); int n = argc - 3;
    ll *v = xmalloc(sizeof(ll) * n); for (int i = 0; i < n; i++) v[i] = atoll(argv[3 + i]);
    int *a = xmalloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) { int r = 0; for (int j = 0; j < n; j++) r += v[j] < v[i]; a[i] = r; }
    static const char *nm[NOPS] = {"sa", "sb", "ss", "pa", "pb", "ra", "rb", "rr", "rra", "rrb", "rrr"};
    OpVec ops; ov_init(&ops); char w[8];
    while (fscanf(f, "%7s", w) == 1) for (int o = 0; o < NOPS; o++) if (!strcmp(w, nm[o])) ov_push(&ops, o);
    int before = ops.n, tot = 0, s;
    do { s = peep_optimize(a, n, &ops, K); tot += s; } while (s > 0);
    printf("%d %d %d\n", before, ops.n, verify_ops(a, n, &ops));
    return 0;
}
