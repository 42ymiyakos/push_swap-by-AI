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

typedef struct { int nlis, ktop, bw, bws, bc, kper, bwm, nm; int mlist4[5]; int verbose, kdir; } InitCfg;
static const int invop[NOPS] = {SA, SB, SS, PB, PA, RRA, RRB, RRR, RA, RB, RR};
/* ndir = 2: also solve the inverse permutation; its reversed, inverted op sequence sorts a0 */
static void build_init(const int *a0, int n, const InitCfg *C, int ndir, OpVec *base, int *lis_out) {
    int nlis = C->nlis, ktop = C->ktop, bw = C->bw, bws = C->bws, bc = C->bc, kper = C->kper, bwm = C->bwm, nm = C->nm;
    const int *mlist4 = C->mlist4;
    int *arr[2]; arr[0] = (int *)a0; arr[1] = NULL;
    if (ndir == 2) { arr[1] = xmalloc(sizeof(int) * n); for (int i = 0; i < n; i++) arr[1][a0[i]] = i; }
    int starts = imin(n, nlis), nK[2] = {0, 0};
    char **Ks[2] = {NULL, NULL};
    SC *st2 = xmalloc(sizeof(SC) * (2 * ndir * nm * kper + 1));
    int *st2m = xmalloc(sizeof(int) * (2 * ndir * nm * kper + 1)), n2 = 0;
    OpVec tmp; ov_init(&tmp);
    double tb0 = now_sec();
    for (int dir = 0; dir < ndir; dir++) {
        const int *a = arr[dir];
        /* candidate K sets: LIS from many starts (distinct), screened with a small beam */
        Ks[dir] = xmalloc(sizeof(char *) * starts);
        for (int t = 0; t < starts; t++) {
            int st = (int)((long long)t * n / starts);
            char *K = xmalloc(n);
            lis_from(a, n, st, K);
            int dup = 0;
            for (int k = 0; k < nK[dir]; k++) if (!memcmp(Ks[dir][k], K, n)) { dup = 1; break; }
            if (dup) free(K); else Ks[dir][nK[dir]++] = K;
        }
        /* stage 1: greedy screening per chunk count; stage 2: small beam (kper per count); stage 3: full beam */
        SC *sc = xmalloc(sizeof(SC) * nK[dir]);
        for (int cm = 0; cm < 2; cm++) for (int mi = 0; mi < nm; mi++) {
            int mm = mlist4[mi]; g_cmode = cm;
            for (int k = 0; k < nK[dir]; k++) { chunk_construct(a, n, Ks[dir][k], mm, bws, bc, &tmp); sc[k].size = tmp.n; sc[k].k = k; }
            qsort(sc, nK[dir], sizeof(SC), cmp_sc);
            for (int k = 0; k < imin(nK[dir], kper); k++) {
                int sz = sc[k].size;
                if (bwm > 0) { chunk_construct(a, n, Ks[dir][sc[k].k], mm, bwm, bc, &tmp); sz = tmp.n; }
                st2[n2].size = sz; st2[n2].k = sc[k].k; st2m[n2] = mm + 1000 * cm + 10000 * dir; n2++;
            }
        }
        free(sc);
        if (C->verbose) fprintf(stderr, "dir %d screened t=%.2f\n", dir, now_sec() - tb0);
    }
    /* sort (size, k, m) */
    for (int x = 1; x < n2; x++) {
        SC e = st2[x]; int em = st2m[x], y = x - 1;
        while (y >= 0 && (st2[y].size > e.size || (st2[y].size == e.size && (st2[y].k > e.k || (st2[y].k == e.k && st2m[y] > em))))) {
            st2[y + 1] = st2[y]; st2m[y + 1] = st2m[y]; y--;
        }
        st2[y + 1] = e; st2m[y + 1] = em;
    }
    int haveBase = 0, bestDir = 0, bestK = 0;
    int taken[2] = {0, 0}, done = 0;
    for (int k = 0; k < n2 && done < ktop; k++) {
        int dir = st2m[k] / 10000;
        if (C->kdir > 0) { if (taken[dir] >= C->kdir) continue; taken[dir]++; done = (taken[0] >= C->kdir) + (taken[1] >= C->kdir) >= ndir ? ktop : 0; }
        else done++;
        g_cmode = (st2m[k] / 1000) % 10;
        chunk_construct(arr[dir], n, Ks[dir][st2[k].k], st2m[k] % 1000, bw, bc, &tmp);
        if (dir) {
            OpVec r; ov_init(&r);
            for (int i = tmp.n - 1; i >= 0; i--) ov_push(&r, invop[tmp.d[i]]);
            ov_copy(&tmp, &r); ov_free(&r);
        }
        if (C->verbose) fprintf(stderr, "full dir=%d stage2=%d full=%d t=%.2f\n", dir, st2[k].size, tmp.n, now_sec() - tb0);
        if (!haveBase || tmp.n < base->n) { ov_copy(base, &tmp); haveBase = 1; bestDir = dir; bestK = st2[k].k; }
    }
    if (C->verbose) fprintf(stderr, "init dir=%d\n", bestDir);
    free(st2); free(st2m);
    int lis_len = 0;
    for (int i = 0; i < n; i++) lis_len += Ks[bestDir][bestK][i];
    for (int dir = 0; dir < ndir; dir++) { for (int k = 0; k < nK[dir]; k++) free(Ks[dir][k]); free(Ks[dir]); }
    ov_free(&tmp);
    if (arr[1]) free(arr[1]);
    *lis_out = lis_len;
}

int main(int argc, char **argv) {
    Params P = {-1, -1, 1, 0, 3, 8, 6, 6, 1, 0, 0.5, 1, 1, 1, 16, 0, 1, 200, 2, 0, -1};
    int first = 1, seglenSet = 0, ktopOpt = -1, nlisOpt = -1, kdirOpt = -1, dualig = 1, bwOpt = -1, bwmOpt = -1, bwsOpt = -1;
    for (; first < argc && argv[first][0] == '-' && argv[first][1] == '-'; first++) {
        const char *a = argv[first];
        if (!strncmp(a, "--time=", 7)) P.time = atof(a + 7);
        else if (!strncmp(a, "--iters=", 8)) P.iters = atoll(a + 8);
        else if (!strncmp(a, "--seed=", 7)) P.seed = strtoull(a + 7, NULL, 10);
        else if (!strncmp(a, "--verbose=", 10)) P.verbose = atoi(a + 10);
        else if (!strncmp(a, "--a2a=", 6)) P.a2a = atoi(a + 6);
        else if (!strncmp(a, "--ubw=", 6)) P.ubw = atoi(a + 6);
        else if (!strncmp(a, "--irot=", 7)) P.irot = atoi(a + 7);
        else if (!strncmp(a, "--dbo=", 6)) P.dbo = atoi(a + 6);
        else if (!strncmp(a, "--dbias=", 8)) P.dbias = atoi(a + 8);
        else if (!strncmp(a, "--kick=", 7)) P.kick = atoi(a + 7);
        else if (!strncmp(a, "--kickmul=", 10)) P.kickmul = atoi(a + 10);
        else if (!strncmp(a, "--kickbest=", 11)) P.kickbest = atoi(a + 11);
        else if (!strncmp(a, "--nlis=", 7)) nlisOpt = atoi(a + 7);
        else if (!strncmp(a, "--dualig=", 9)) dualig = atoi(a + 9);
        else if (!strncmp(a, "--kdir=", 7)) kdirOpt = atoi(a + 7);
        else if (!strncmp(a, "--ktop=", 7)) ktopOpt = atoi(a + 7);
        else if (!strncmp(a, "--bw=", 5)) bwOpt = atoi(a + 5);
        else if (!strncmp(a, "--bwm=", 6)) bwmOpt = atoi(a + 6);
        else if (!strncmp(a, "--bws=", 6)) bwsOpt = atoi(a + 6);
        else if (!strncmp(a, "--pqd=", 6)) P.pqd = atof(a + 6);
        else if (!strncmp(a, "--seglen=", 9)) { P.seglen = atoi(a + 9); seglenSet = 1; }
        else if (!strncmp(a, "--dual=", 7)) P.dual = atoi(a + 7);
        else if (!strncmp(a, "--a2amax=", 9)) P.a2amax = atoi(a + 9);
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
        g_cmode = 0;
        chunk_construct(a0, n, K, m, 0, 2, &ops);
        print_ops(&ops); ov_free(&ops); free(K); free(a0); return 0;
    }
    if (P.time < 0) P.time = n <= 150 ? 15.0 : (n <= 300 ? 25.0 : 45.0);
    int nlis = 50, ktop = 3, bw = 1000, bws = 20, bc = 2, kper = 3, bwm = 100;
    int mlist4[5] = {2, 3, 4, 0, 0}, nm = 3;
    if (n > 150 && n <= 300) { int t[4] = {3, 4, 5, 6}; memcpy(mlist4, t, sizeof t); nm = 4; }
    if (n > 300) { if (!seglenSet) P.seglen = 16; P.rmax = 9; P.vwin = 9; nlis = 30; ktop = 3; bw = 500; bws = 0; bwm = 20; int t[5] = {4, 6, 8, 10, 12}; memcpy(mlist4, t, sizeof t); nm = 5; }
    if (bwOpt > 0) bw = bwOpt;
    if (ktopOpt > 0) ktop = ktopOpt;
    if (nlisOpt > 0) nlis = nlisOpt;
    if (bwmOpt > 0) bwm = bwmOpt;
    if (bwsOpt >= 0) bws = bwsOpt;
    int peepK[3] = {4, 6, 8}, npeep = n <= 150 ? 3 : 2;
    double t0g = now_sec();

    InitCfg C = {nlis, ktop, bw, bws, bc, kper, bwm, nm, {0}, 0, 0};
    memcpy(C.mlist4, mlist4, sizeof C.mlist4);
    OpVec base; ov_init(&base);
    int lis_len = 0;
    /* dual (inverse-permutation) initial solutions, ktop full beams per direction: on by default */
    if (P.dual < 0) P.dual = 1;
    if (kdirOpt < 0) kdirOpt = P.dual ? ktop : 0;
    C.verbose = P.verbose; C.kdir = kdirOpt;
    build_init(a0, n, &C, P.dual ? 2 : 1, &base, &lis_len);
    Rep rep; memset(&rep, 0, sizeof rep);
    if (!verify_ops(a0, n, &base) || !convert(a0, n, &base, &rep)) {
        /* should not happen: fall back to the constructor output */
        print_ops(&base); return 0;
    }
    if (P.verbose) fprintf(stderr, "init=%d t=%.2f\n", base.n, now_sec() - t0g);
    Rng rng; rng_seed(&rng, P.seed + 99);
    double t0 = now_sec();
    P.time = P.time - (t0 - t0g) - (n <= 150 ? 0.8 : 1.0);  /* reserve for peephole (measured ~0.5 s) */
    if (P.time < 0.5) P.time = 0.5;
    Rep best; memset(&best, 0, sizeof best);
    int bestCost; long long iters;
    if (dualig <= 1) {
        ig_run(&rep, &P, &rng, &best, &bestCost, &iters);
        emit(&best, &ops);
        if (ops.n != bestCost || !verify_ops(a0, n, &ops)) {
            fprintf(stderr, "internal: emit mismatch %d vs %d\n", ops.n, bestCost);
            ov_copy(&ops, &base);
        }
    } else {
        /* alternate IG phases between P and its inverse (reversed, inverted op sequences) */
        int *inv = xmalloc(sizeof(int) * n);
        for (int i = 0; i < n; i++) inv[a0[i]] = i;
        int dir = 0; iters = 0;
        double T = P.time, tst = now_sec();
        OpVec cur; ov_init(&cur);
        g_convSA = 1;
        for (int ph = 0; ph < dualig; ph++) {
            Params Pp = P; Pp.time = T * (ph + 1) / dualig - (now_sec() - tst);
            if (Pp.time < 0.05) Pp.time = 0.05;
            long long it = 0;
            ig_run(&rep, &Pp, &rng, &best, &bestCost, &it); iters += it;
            emit(&best, &cur);
            if (P.verbose) fprintf(stderr, "phase %d dir %d cost %d iters %lld\n", ph, dir, cur.n, it);
            if (ph + 1 == dualig) break;
            OpVec r2; ov_init(&r2);
            for (int i = cur.n - 1; i >= 0; i--) ov_push(&r2, invop[cur.d[i]]);
            Rep nr; memset(&nr, 0, sizeof nr);
            const int *other = dir ? a0 : inv;
            if (verify_ops(other, n, &r2) && convert(other, n, &r2, &nr)) {
                int fc = full_cost(&nr);
                if (P.verbose) fprintf(stderr, "  convert ok: %d -> %d\n", cur.n, fc);
                rep_free(&rep); rep = nr; dir ^= 1;
            } else {
                if (P.verbose) fprintf(stderr, "  convert failed\n");
                rep_free(&nr); rep_copy(&rep, &best);
            }
            ov_free(&r2);
        }
        g_convSA = 0;
        if (dir) { ov_free(&ops); ov_init(&ops); for (int i = cur.n - 1; i >= 0; i--) ov_push(&ops, invop[cur.d[i]]); }
        else ov_copy(&ops, &cur);
        if (!verify_ops(a0, n, &ops)) { fprintf(stderr, "internal: dualig verify failed\n"); ov_copy(&ops, &base); }
        ov_free(&cur); free(inv);
    }
    int saved = 0;
    double tp0 = now_sec();
    for (int k = 0; k < npeep; k++) saved += peep_optimize(a0, n, &ops, peepK[k]);
    double tpeep = now_sec() - tp0;
    if (!verify_ops(a0, n, &ops)) emit(&best, &ops);
    if (P.verbose)
        fprintf(stderr, "n=%d lis=%d constructor=%d final=%d peep=%d (%.2fs) time=%.2f iters=%lld\n", n, lis_len, base.n, ops.n, saved, tpeep,
                now_sec() - t0, iters);
    print_ops(&ops);
    ov_free(&ops); ov_free(&base); rep_free(&rep); rep_free(&best); free(a0);
    return 0;
}
