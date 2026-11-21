/* push_swap (C port of PSIG). Usage: push_swap [--time=S] [--iters=N] [--seed=S] [--verbose=1] numbers... */
#include "ps.h"

static int parse_args(int argc, char **argv, int first, ll **out, int *nout) {
    int cap = 16, n = 0;
    ll *v = xmalloc(sizeof(ll) * cap);
    for (int i = first; i < argc; i++) {
        const char *s = argv[i];
        while (*s) {
            while (*s == ' ' || *s == '\t') s++;
            if (!*s) break;
            const char *st = s; char *end;
            if (*s == '+' || *s == '-') s++;
            if (!(*s >= '0' && *s <= '9')) { free(v); return 0; }
            errno = 0;
            long long x = strtoll(st, &end, 10);
            if (errno || x < INT_MIN || x > INT_MAX) { free(v); return 0; }
            if (*end && *end != ' ' && *end != '\t') { free(v); return 0; }
            if (n == cap) { cap *= 2; v = xrealloc(v, sizeof(ll) * cap); }
            v[n++] = x; s = end;
        }
    }
    *out = v; *nout = n;
    return 1;
}
static int cmp_ll(const void *a, const void *b) { ll x = *(const ll *)a, y = *(const ll *)b; return (x > y) - (x < y); }

typedef struct { int size, k; } SC;
static int cmp_sc(const void *a, const void *b) {
    const SC *x = a, *y = b;
    if (x->size != y->size) return x->size < y->size ? -1 : 1;
    return (x->k > y->k) - (x->k < y->k);
}

int main(int argc, char **argv) {
    Params P = {-1, -1, 1, 0, 3, 8, 6, 6, 1, 0, 0.5};
    int first = 1;
    for (; first < argc && argv[first][0] == '-' && argv[first][1] == '-'; first++) {
        const char *a = argv[first];
        if (!strncmp(a, "--time=", 7)) P.time = atof(a + 7);
        else if (!strncmp(a, "--iters=", 8)) P.iters = atoll(a + 8);
        else if (!strncmp(a, "--seed=", 7)) P.seed = strtoull(a + 7, NULL, 10);
        else if (!strncmp(a, "--verbose=", 10)) P.verbose = atoi(a + 10);
        else { fprintf(stderr, "Error\n"); return 1; }
    }
    ll *vals; int n;
    if (!parse_args(argc, argv, first, &vals, &n)) { fprintf(stderr, "Error\n"); return 1; }
    ll *sv = xmalloc(sizeof(ll) * (n + 1));
    memcpy(sv, vals, sizeof(ll) * n);
    qsort(sv, n, sizeof(ll), cmp_ll);
    for (int i = 1; i < n; i++) if (sv[i] == sv[i - 1]) { fprintf(stderr, "Error\n"); return 1; }
    if (n == 0) return 0;
    int *a0 = xmalloc(sizeof(int) * n);
    int sorted = 1;
    for (int i = 0; i < n; i++) {
        int lo = 0, hi = n;
        while (lo < hi) { int mid = (lo + hi) >> 1; if (sv[mid] < vals[i]) lo = mid + 1; else hi = mid; }
        a0[i] = lo;
        if (a0[i] != i) sorted = 0;
    }
    free(sv); free(vals);
    if (sorted) { free(a0); return 0; }
    OpVec ops; ov_init(&ops);
    if (n <= 8) { exact_solve(a0, n, &ops); print_ops(&ops); ov_free(&ops); free(a0); return 0; }
    if (n > 1000) {  /* values are encoded in 10 bits: constructive fallback */
        Rng rr; rng_seed(&rr, 1);
        char *K = xmalloc(n);
        circular_lis(a0, n, &rr, K);
        int m = imax(2, n / 80); if (m > 120) m = 120;
        chunk_construct(a0, n, K, m, 0, 2, &ops);
        print_ops(&ops); ov_free(&ops); free(K); free(a0); return 0;
    }
    if (P.time < 0) P.time = n <= 150 ? 10.0 : (n <= 300 ? 20.0 : 40.0);
    int nlis = 50, ktop = 3, bw = 1000, bws = 20, bc = 2;
    if (n > 300) { P.seglen = 16; P.rmax = 9; P.vwin = 9; nlis = 30; ktop = 3; bw = 100; bws = 0; }
    int peepK[3] = {4, 6, 8}, npeep = n <= 150 ? 3 : 2;
    double t0g = now_sec();
    int m = n <= 150 ? 2 : (n <= 300 ? 4 : 6);

    /* candidate K sets: LIS from many starts (distinct), screened with a small beam */
    int starts = imin(n, nlis), nK = 0;
    char **Ks = xmalloc(sizeof(char *) * starts);
    for (int t = 0; t < starts; t++) {
        int st = (int)((long long)t * n / starts);
        char *K = xmalloc(n);
        lis_from(a0, n, st, K);
        int dup = 0;
        for (int k = 0; k < nK; k++) if (!memcmp(Ks[k], K, n)) { dup = 1; break; }
        if (dup) free(K); else Ks[nK++] = K;
    }
    SC *sc = xmalloc(sizeof(SC) * nK);
    OpVec tmp; ov_init(&tmp);
    for (int k = 0; k < nK; k++) { chunk_construct(a0, n, Ks[k], m, bws, bc, &tmp); sc[k].size = tmp.n; sc[k].k = k; }
    qsort(sc, nK, sizeof(SC), cmp_sc);
    OpVec base; ov_init(&base);
    int haveBase = 0;
    for (int k = 0; k < imin(nK, ktop); k++) {
        chunk_construct(a0, n, Ks[sc[k].k], m, bw, bc, &tmp);
        if (!haveBase || tmp.n < base.n) { ov_copy(&base, &tmp); haveBase = 1; }
    }
    int lis_len = 0;
    for (int i = 0; i < n; i++) lis_len += Ks[sc[0].k][i];
    for (int k = 0; k < nK; k++) free(Ks[k]);
    free(Ks); free(sc); ov_free(&tmp);

    Rep rep; memset(&rep, 0, sizeof rep);
    if (!verify_ops(a0, n, &base) || !convert(a0, n, &base, &rep)) {
        /* should not happen: fall back to the constructor output */
        print_ops(&base); return 0;
    }
    if (P.verbose) fprintf(stderr, "init=%d t=%.2f\n", base.n, now_sec() - t0g);
    Rng rng; rng_seed(&rng, P.seed + 99);
    double t0 = now_sec();
    P.time = P.time - (t0 - t0g) - (n <= 150 ? 2.0 : 4.0);
    if (P.time < 0.5) P.time = 0.5;
    Rep best; memset(&best, 0, sizeof best);
    int bestCost; long long iters;
    ig_run(&rep, &P, &rng, &best, &bestCost, &iters);
    emit(&best, &ops);
    if (ops.n != bestCost || !verify_ops(a0, n, &ops)) {
        fprintf(stderr, "internal: emit mismatch %d vs %d\n", ops.n, bestCost);
        ov_copy(&ops, &base);
    }
    int saved = 0;
    for (int k = 0; k < npeep; k++) saved += peep_optimize(a0, n, &ops, peepK[k]);
    if (!verify_ops(a0, n, &ops)) emit(&best, &ops);
    if (P.verbose)
        fprintf(stderr, "n=%d lis=%d constructor=%d final=%d peep=%d time=%.2f iters=%lld\n", n, lis_len, base.n, ops.n, saved,
                now_sec() - t0, iters);
    print_ops(&ops);
    ov_free(&ops); ov_free(&base); rep_free(&rep); rep_free(&best); free(a0);
    return 0;
}
