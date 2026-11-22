/* push_swap (C port of PSIG). Usage: push_swap [--time=S] [--iters=N] [--seed=S] [--verbose=1] numbers... */
#include "ps.h"

static int g_convSlack = 1000000;
/* a converted representation is used only if it emits a valid sequence of the computed cost */
static int conv_ok(const int *a, int n, const Rep *r, int fc) {
    if (!g_convFix) return 1;
    OpVec o; ov_init(&o); emit(r, &o);
    int ok = o.n == fc && verify_ops(a, n, &o);
    ov_free(&o); return ok;
}
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

typedef struct { int nlis, ktop, bw, bws, bc, kper, bwm, nm; int mlist4[5]; int verbose, kdir, isswap, swratio, force; } InitCfg;
#define MAXCAND 16
static OpVec g_cand[MAXCAND]; static int g_ncand = 0;  /* all full-beam start solutions (for probing) */
static const int invop[NOPS] = {SA, SB, SS, PB, PA, RRA, RRB, RRR, RA, RB, RR};
/* ndir = 2: also solve the inverse permutation; its reversed, inverted op sequence sorts a0 */
/* ndir = 2: also solve the inverse permutation; its reversed, inverted op sequence sorts a0.
 * Candidate groups ("virtual directions"): plain LIS for each direction, and with isswap, IS+swap candidates
 * for each direction whose best IS+swap length is >= swratio/10 x its LIS length (structured inputs such as
 * descending blocks; random inputs never pass). Each group is screened separately and gets kdir full beams. */
static void build_init(const int *a0, int n, const InitCfg *C, int ndir, OpVec *base, int *lis_out) {
    int nlis = C->nlis, ktop = C->ktop, bw = C->bw, bws = C->bws, bc = C->bc, kper = C->kper, bwm = C->bwm, nm = C->nm;
    const int *mlist4 = C->mlist4;
    int *arr[2]; arr[0] = (int *)a0; arr[1] = NULL;
    if (ndir == 2) { arr[1] = xmalloc(sizeof(int) * n); for (int i = 0; i < n; i++) arr[1][a0[i]] = i; }
    int nv = C->isswap ? 2 * ndir : ndir;
    int starts = imin(n, nlis), nK[4] = {0, 0, 0, 0}, maxLis[2] = {0, 0};
    char **Ks[4] = {NULL, NULL, NULL, NULL};
    int **Pr[4] = {NULL, NULL, NULL, NULL};  /* swap pairs per candidate (NULL = plain LIS) */
    SC *st2 = xmalloc(sizeof(SC) * (2 * nv * nm * kper + 1));
    int *st2m = xmalloc(sizeof(int) * (2 * nv * nm * kper + 1)), n2 = 0;
    OpVec tmp; ov_init(&tmp);
    double tb0 = now_sec();
    for (int vd = 0; vd < nv; vd++) {
        int dir = vd % ndir, sw = vd >= ndir;
        const int *a = arr[dir];
        Ks[vd] = xmalloc(sizeof(char *) * starts); Pr[vd] = xcalloc(starts, sizeof(int *));
        int maxSw = 0;
        for (int t = 0; t < starts; t++) {
            int st = (int)((long long)t * n / starts);
            char *K = xmalloc(n); int *pr = NULL;
            if (!sw) { lis_from(a, n, st, K); int L = 0; for (int v = 0; v < n; v++) L += K[v]; if (L > maxLis[dir]) maxLis[dir] = L; }
            else {
                pr = xmalloc(sizeof(int) * n);
                int L = isswap_from(a, n, st, K, pr);
                int any = 0; for (int v = 0; v < n; v++) if (pr[v] >= 0) { any = 1; break; }
                if (!any) { free(K); free(pr); continue; }
                if (L > maxSw) maxSw = L;
            }
            int dup = 0;
            for (int k = 0; k < nK[vd]; k++) if (!memcmp(Ks[vd][k], K, n) && (!pr || !memcmp(pr, Pr[vd][k], sizeof(int) * n))) { dup = 1; break; }
            if (dup) { free(K); free(pr); } else { Pr[vd][nK[vd]] = pr; Ks[vd][nK[vd]++] = K; }
        }
        if (sw && maxSw * 10 < maxLis[dir] * C->swratio) {  /* group not used for this input */
            for (int k = 0; k < nK[vd]; k++) { free(Ks[vd][k]); free(Pr[vd][k]); }
            nK[vd] = 0;
        }
        if (!nK[vd]) continue;
        /* stage 1: greedy screening per chunk count; stage 2: small beam (kper per count); stage 3: full beam */
        SC *sc = xmalloc(sizeof(SC) * nK[vd]);
        for (int cm = 0; cm < 2; cm++) for (int mi = 0; mi < nm; mi++) {
            int mm = mlist4[mi]; g_cmode = cm;
            for (int k = 0; k < nK[vd]; k++) { g_pair = Pr[vd][k]; chunk_construct(a, n, Ks[vd][k], mm, bws, bc, &tmp); sc[k].size = tmp.n; sc[k].k = k; }
            qsort(sc, nK[vd], sizeof(SC), cmp_sc);
            for (int k = 0; k < imin(nK[vd], kper); k++) {
                int sz = sc[k].size;
                if (bwm > 0) { g_pair = Pr[vd][sc[k].k]; chunk_construct(a, n, Ks[vd][sc[k].k], mm, bwm, bc, &tmp); sz = tmp.n; }
                st2[n2].size = sz; st2[n2].k = sc[k].k; st2m[n2] = mm + 1000 * cm + 10000 * vd; n2++;
            }
        }
        free(sc);
        if (C->verbose) fprintf(stderr, "group %d (%d candidates) screened t=%.2f\n", vd, nK[vd], now_sec() - tb0);
    }
    /* sort (size, k, m) */
    for (int x = 1; x < n2; x++) {
        SC e = st2[x]; int em = st2m[x], y = x - 1;
        while (y >= 0 && (st2[y].size > e.size || (st2[y].size == e.size && (st2[y].k > e.k || (st2[y].k == e.k && st2m[y] > em))))) {
            st2[y + 1] = st2[y]; st2m[y + 1] = st2m[y]; y--;
        }
        st2[y + 1] = e; st2m[y + 1] = em;
    }
    int haveBase = 0, bestV = 0, bestK = 0;
    int taken[4] = {0, 0, 0, 0}, done = 0, ngroups = 0, nfull = 0;
    for (int vd = 0; vd < nv; vd++) ngroups += nK[vd] > 0;
    for (int k = 0; k < n2 && done < ktop; k++) {
        int vd = st2m[k] / 10000, dir = vd % ndir;
        if (C->kdir > 0) {
            if (taken[vd] >= C->kdir) continue;
            taken[vd]++;
            int full = 0; for (int g = 0; g < nv; g++) full += nK[g] > 0 && taken[g] >= C->kdir;
            done = full >= ngroups ? ktop : 0;
        } else done++;
        g_cmode = (st2m[k] / 1000) % 10;
        g_pair = Pr[vd][st2[k].k];
        chunk_construct(arr[dir], n, Ks[vd][st2[k].k], st2m[k] % 1000, bw, bc, &tmp);
        if (dir) {
            OpVec r; ov_init(&r);
            for (int i = tmp.n - 1; i >= 0; i--) ov_push(&r, invop[tmp.d[i]]);
            ov_copy(&tmp, &r); ov_free(&r);
        }
        if (C->verbose) fprintf(stderr, "full group=%d stage2=%d full=%d t=%.2f\n", vd, st2[k].size, tmp.n, now_sec() - tb0);
        if (g_ncand < MAXCAND) { ov_init(&g_cand[g_ncand]); ov_copy(&g_cand[g_ncand], &tmp); g_ncand++; }
        if (C->force >= 0 ? nfull == C->force : (!haveBase || tmp.n < base->n)) { ov_copy(base, &tmp); haveBase = 1; bestV = vd; bestK = st2[k].k; }
        if (!haveBase) { ov_copy(base, &tmp); haveBase = 1; bestV = vd; bestK = st2[k].k; }
        nfull++;
    }
    if (C->verbose) fprintf(stderr, "init group=%d\n", bestV);
    g_pair = NULL;
    free(st2); free(st2m);
    int lis_len = 0;
    for (int i = 0; i < n; i++) lis_len += Ks[bestV][bestK][i];
    for (int vd = 0; vd < nv; vd++) { for (int k = 0; k < nK[vd]; k++) { free(Ks[vd][k]); free(Pr[vd][k]); } free(Ks[vd]); free(Pr[vd]); }
    ov_free(&tmp);
    if (arr[1]) free(arr[1]);
    *lis_out = lis_len;
}

int main(int argc, char **argv) {
    Params P = {-1, -1, 1, 0, 3, 8, 6, 6, 1, 0, 0.5, 1, 1, 1, 16, 0, 1, -1, 2, 0, -1};
    int first = 1, seglenSet = 0, ktopOpt = -1, nlisOpt = -1, kdirOpt = -1, dualig = 8, isswapOpt = 1, swrOpt = 0, forceOpt = -1, probe = 0, bwOpt = -1, bwmOpt = -1, bwsOpt = -1;
    double probet = 2.0;
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
        else if (!strncmp(a, "--lahc=", 7)) P.lahc = atoi(a + 7);
        else if (!strncmp(a, "--kickmul=", 10)) P.kickmul = atoi(a + 10);
        else if (!strncmp(a, "--kickbest=", 11)) P.kickbest = atoi(a + 11);
        else if (!strncmp(a, "--nlis=", 7)) nlisOpt = atoi(a + 7);
        else if (!strncmp(a, "--dualig=", 9)) dualig = atoi(a + 9);
        else if (!strncmp(a, "--probe=", 8)) probe = atoi(a + 8);
        else if (!strncmp(a, "--probet=", 9)) probet = atof(a + 9);
        else if (!strncmp(a, "--forceinit=", 12)) forceOpt = atoi(a + 12);
        else if (!strncmp(a, "--swratio=", 10)) swrOpt = atoi(a + 10);
        else if (!strncmp(a, "--isswap=", 9)) isswapOpt = atoi(a + 9);
        else if (!strncmp(a, "--kdir=", 7)) kdirOpt = atoi(a + 7);
        else if (!strncmp(a, "--ktop=", 7)) ktopOpt = atoi(a + 7);
        else if (!strncmp(a, "--bw=", 5)) bwOpt = atoi(a + 5);
        else if (!strncmp(a, "--bwm=", 6)) bwmOpt = atoi(a + 6);
        else if (!strncmp(a, "--bws=", 6)) bwsOpt = atoi(a + 6);
        else if (!strncmp(a, "--pqd=", 6)) P.pqd = atof(a + 6);
        else if (!strncmp(a, "--seglen=", 9)) { P.seglen = atoi(a + 9); seglenSet = 1; }
        else if (!strncmp(a, "--dual=", 7)) P.dual = atoi(a + 7);
        else if (!strncmp(a, "--convfix=", 10)) g_convFix = atoi(a + 10);
        else if (!strncmp(a, "--convslack=", 12)) g_convSlack = atoi(a + 12);
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
    if (P.kick < 0) P.kick = P.time > 90.0 ? 50 : 200;  /* long budgets: kick more often (log 72) */
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

    InitCfg C = {nlis, ktop, bw, bws, bc, kper, bwm, nm, {0}, 0, 0, 0, 16, -1};
    memcpy(C.mlist4, mlist4, sizeof C.mlist4);
    OpVec base; ov_init(&base);
    int lis_len = 0;
    /* dual (inverse-permutation) initial solutions, ktop full beams per direction: on by default */
    if (P.dual < 0) P.dual = 1;
    if (kdirOpt < 0) kdirOpt = P.dual ? ktop : 0;
    C.verbose = P.verbose; C.kdir = kdirOpt; C.isswap = isswapOpt; C.force = forceOpt; if (swrOpt > 0) C.swratio = swrOpt;
    build_init(a0, n, &C, P.dual ? 2 : 1, &base, &lis_len);
    Rep rep; memset(&rep, 0, sizeof rep);
    g_convSA = 1;  /* the base may contain SA (IS+swap starts) */
    int convOK = verify_ops(a0, n, &base) && convert(a0, n, &base, &rep);
    g_convSA = 0;
    if (!convOK) {
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
    if (probe > 1 && g_ncand > 1) {
        /* short IG probe on the best few start solutions; continue from the probe that got furthest */
        int ord[MAXCAND];
        for (int i = 0; i < g_ncand; i++) ord[i] = i;
        for (int i = 1; i < g_ncand; i++) { int e = ord[i], j = i - 1; while (j >= 0 && g_cand[ord[j]].n > g_cand[e].n) { ord[j + 1] = ord[j]; j--; } ord[j + 1] = e; }
        int np = imin(probe, g_ncand), bestP = INT_MAX;
        double tps = now_sec();
        Rep pr; memset(&pr, 0, sizeof pr);
        for (int q = 0; q < np; q++) {
            Rep cr; memset(&cr, 0, sizeof cr);
            g_convSA = 1;
            int ok = convert(a0, n, &g_cand[ord[q]], &cr);
            g_convSA = 0;
            if (!ok) { rep_free(&cr); continue; }
            Params Pq = P; Pq.time = probet;
            int c; long long it;
            ig_run(&cr, &Pq, &rng, &pr, &c, &it);
            if (P.verbose) fprintf(stderr, "probe %d: start %d -> %d\n", q, g_cand[ord[q]].n, c);
            if (c < bestP) { bestP = c; rep_copy(&rep, &pr); }
            rep_free(&cr);
        }
        rep_free(&pr);
        P.time -= now_sec() - tps;
        if (P.time < 0.5) P.time = 0.5;
    }
    for (int i = 0; i < g_ncand; i++) ov_free(&g_cand[i]);
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
        OpVec gbest; ov_init(&gbest);
        g_convSA = 1;
        for (int ph = 0; ph < dualig; ph++) {
            Params Pp = P; Pp.time = T * (ph + 1) / dualig - (now_sec() - tst);
            if (Pp.time < 0.05) Pp.time = 0.05;
            long long it = 0;
            ig_run(&rep, &Pp, &rng, &best, &bestCost, &it); iters += it;
            emit(&best, &cur);
            if (P.verbose) fprintf(stderr, "phase %d dir %d cost %d iters %lld\n", ph, dir, cur.n, it);
            if (gbest.n == 0 || cur.n < gbest.n) {
                ov_free(&gbest); ov_init(&gbest);
                if (dir) for (int i = cur.n - 1; i >= 0; i--) ov_push(&gbest, invop[cur.d[i]]);
                else ov_copy(&gbest, &cur);
            }
            if (ph + 1 == dualig) break;
            OpVec r2; ov_init(&r2);
            for (int i = cur.n - 1; i >= 0; i--) ov_push(&r2, invop[cur.d[i]]);
            Rep nr; memset(&nr, 0, sizeof nr);
            const int *other = dir ? a0 : inv;
            int fc = -1;
            if (verify_ops(other, n, &r2) && convert(other, n, &r2, &nr) && (fc = full_cost(&nr)) <= cur.n + g_convSlack && conv_ok(other, n, &nr, fc)) {
                if (P.verbose) fprintf(stderr, "  convert ok: %d -> %d (fixed %d)\n", cur.n, fc, g_convFixed);
                rep_free(&rep); rep = nr; dir ^= 1;
            } else {
                if (P.verbose) fprintf(stderr, "  convert failed\n");
                rep_free(&nr); rep_copy(&rep, &best);
            }
            ov_free(&r2);
        }
        g_convSA = 0;
        ov_copy(&ops, &gbest); ov_free(&gbest);
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
