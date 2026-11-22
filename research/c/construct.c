/* Initial solutions: LIS choice, chunk A->B (middle-out top/bottom chunks), B->A greedy or beam. */
#include "ps.h"

/* LIS reading A cyclically from position s; mask[value] = 1 for members. */
static int lis_core(const int *a, int n, int s, int *prev, int *tails, int *tailidx) {
    int L = 0;
    for (int k = 0; k < n; k++) {
        int i = (s + k) % n, x = a[i];
        int p = lb_int(tails, L, x);
        prev[i] = p ? tailidx[p - 1] : -1;
        tails[p] = x; tailidx[p] = i;
        if (p == L) L++;
    }
    return L;
}
void lis_from(const int *a, int n, int s, char *mask) {
    int *prev = xmalloc(sizeof(int) * n), *tails = xmalloc(sizeof(int) * n), *ti = xmalloc(sizeof(int) * n);
    int L = lis_core(a, n, s, prev, tails, ti);
    memset(mask, 0, n);
    for (int i = ti[L - 1]; i >= 0; i = prev[i]) mask[a[i]] = 1;
    free(prev); free(tails); free(ti);
}
void circular_lis(const int *a, int n, Rng *rng, char *mask) {
    int *prev = xmalloc(sizeof(int) * n), *tails = xmalloc(sizeof(int) * n), *ti = xmalloc(sizeof(int) * n);
    int bestlen = -1, ties = 0;
    memset(mask, 0, n);
    for (int s = 0; s < n; s++) {
        int L = lis_core(a, n, s, prev, tails, ti), take = 0;
        if (L > bestlen) { bestlen = L; ties = 1; take = 1; }
        else if (L == bestlen) { ties++; take = rng_below(rng, ties) == 0; }
        if (take) { memset(mask, 0, n); for (int i = ti[L - 1]; i >= 0; i = prev[i]) mask[a[i]] = 1; }
    }
    free(prev); free(tails); free(ti);
}


/* IS+swap reading A cyclically from position s: longest subsequence that becomes increasing after swapping
 * some disjoint pairs of consecutive members. pairOf[x] = y for a pair (x before y in A, y < x), -1 otherwise.
 * Returns the length. Port of isswap_from in src/construct.hpp. */
typedef struct { int len, id; } FW;
static int fw_gt(FW a, FW b) { return a.len != b.len ? a.len > b.len : a.id > b.id; }
int isswap_from(const int *a, int n, int s, char *mask, int *pairOf) {
    int *b = xmalloc(sizeof(int) * n);
    for (int k = 0; k < n; k++) b[k] = a[(s + k) % n];
    int cap = n * 8 + 16, ns = 0;
    int *slen = xmalloc(sizeof(int) * cap), *spar = xmalloc(sizeof(int) * cap), *sk = xmalloc(sizeof(int) * cap),
        *sp1 = xmalloc(sizeof(int) * cap), *sp2 = xmalloc(sizeof(int) * cap), *nxt = xmalloc(sizeof(int) * cap), *pv = xmalloc(sizeof(int) * cap);
    int *head = xmalloc(sizeof(int) * n);
    for (int k = 0; k < n; k++) head[k] = -1;
    FW *fw = xmalloc(sizeof(FW) * (n + 1));
    for (int i = 0; i <= n; i++) { fw[i].len = 0; fw[i].id = -1; }
    FW best = {0, -1};
#define GROW() if (ns >= cap) { cap *= 2; slen = xrealloc(slen, sizeof(int) * cap); spar = xrealloc(spar, sizeof(int) * cap); sk = xrealloc(sk, sizeof(int) * cap); \
        sp1 = xrealloc(sp1, sizeof(int) * cap); sp2 = xrealloc(sp2, sizeof(int) * cap); nxt = xrealloc(nxt, sizeof(int) * cap); pv = xrealloc(pv, sizeof(int) * cap); }
    for (int t = 0; t < n; t++) {
        FW q = {0, -1};
        for (int x = b[t]; x > 0; x -= x & -x) if (fw_gt(fw[x], q)) q = fw[x];
        GROW();
        int id = ns++; slen[id] = q.len + 1; spar[id] = q.id; sk[id] = 0; sp1[id] = t; sp2[id] = -1;
        for (int k = t + 1; k < n; k++) if (b[k] < b[t]) {
            FW q2 = {0, -1};
            for (int x = b[k]; x > 0; x -= x & -x) if (fw_gt(fw[x], q2)) q2 = fw[x];
            GROW();
            int id2 = ns++; slen[id2] = q2.len + 2; spar[id2] = q2.id; sk[id2] = 1; sp1[id2] = t; sp2[id2] = k;
            pv[id2] = b[t]; nxt[id2] = head[k]; head[k] = id2;
        }
        FW e = {slen[id], id};
        for (int x = b[t] + 1; x <= n; x += x & -x) if (fw_gt(e, fw[x])) fw[x] = e;
        if (fw_gt(e, best)) best = e;
        /* pending pairs ending at t, in insertion order (the C++ version iterates a vector) */
        int cnt = 0;
        for (int p = head[t]; p >= 0; p = nxt[p]) cnt++;
        int *ord = xmalloc(sizeof(int) * (cnt + 1)), c = cnt;
        for (int p = head[t]; p >= 0; p = nxt[p]) ord[--c] = p;
        for (int z = 0; z < cnt; z++) {
            int p = ord[z]; FW e2 = {slen[p], p};
            for (int x = pv[p] + 1; x <= n; x += x & -x) if (fw_gt(e2, fw[x])) fw[x] = e2;
            if (fw_gt(e2, best)) best = e2;
        }
        free(ord);
    }
#undef GROW
    memset(mask, 0, n);
    for (int v = 0; v < n; v++) pairOf[v] = -1;
    for (int id = best.id; id >= 0; id = spar[id]) {
        if (sk[id] == 0) mask[b[sp1[id]]] = 1;
        else { int x = b[sp1[id]], y = b[sp2[id]]; mask[x] = mask[y] = 1; pairOf[x] = y; }
    }
    int L = best.len;
    free(b); free(slen); free(spar); free(sk); free(sp1); free(sp2); free(nxt); free(pv); free(head); free(fw);
    return L;
}

/* with pairs (g_pair != NULL): whenever a pair (x, y) sits at the top of A, swap it */
const int *g_pair = NULL;
typedef struct { Sim s; OpVec *ops; char *swapped; } Builder;
static void bfix(Builder *b) {
    if (!g_pair) return;
    while (b->s.a.len >= 2) {
        int x = dq_get(&b->s.a, 0), y = dq_get(&b->s.a, 1);
        if (g_pair[x] == y && !b->swapped[x]) { b->swapped[x] = 1; sim_apply(&b->s, SA); ov_push(b->ops, SA); }
        else break;
    }
}
static void bop(Builder *b, int o, int k) { for (int i = 0; i < k; i++) { sim_apply(&b->s, o); ov_push(b->ops, o); bfix(b); } }

int g_cmode = 0;  /* 0: equal value ranges, 1: equal counts of moved values */
void chunk_construct(const int *a0, int n, const char *inK, int m, int bw, int bc, OpVec *out) {
    int *crank = xmalloc(sizeof(int) * n), nmv = 0;
    for (int v = 0; v < n; v++) { crank[v] = nmv; if (!inK[v]) nmv++; }
#define CHUNK(v) (g_cmode == 1 && nmv > 0 ? m - 1 - (int)((long long)crank[v] * m / nmv) : m - 1 - (int)((long long)(v) * m / n))
    Builder b; sim_init(&b.s, a0, n); b.ops = out; out->n = 0; b.swapped = g_pair ? xcalloc(n, 1) : NULL;
    bfix(&b);
    int rt[64][2], nr = 0;  /* rounds: (top chunk, bottom chunk) */
    if (m % 2 == 1) { rt[nr][0] = m / 2; rt[nr][1] = -1; nr++; for (int r = 1; r <= m / 2; r++) { rt[nr][0] = m / 2 - r; rt[nr][1] = m / 2 + r; nr++; } }
    else for (int r = 0; r < m / 2; r++) { rt[nr][0] = m / 2 - 1 - r; rt[nr][1] = m / 2 + r; nr++; }
    int pending = 0;
    for (int rd = 0; rd < nr; rd++) {
        int tc = rt[rd][0], bcn = rt[rd][1];
        for (;;) {
            int sa = b.s.a.len, bi = -1, bd = INT_MAX, bdir = 0;
            for (int i = 0; i < sa; i++) {
                int v = dq_get(&b.s.a, i);
                if (inK[v]) continue;
                int c = CHUNK(v);
                if (c != tc && c != bcn) continue;
                int top = c == tc;
                int f = i, cf = top ? imax(f, pending) : f;
                int rv = (sa - i) % sa, cr = rv + (top ? pending : 0);
                int d = imin(cf, cr);
                if (d < bd) { bd = d; bi = i; bdir = cf <= cr ? 0 : 1; }
            }
            if (bi < 0) break;
            int v = dq_get(&b.s.a, bi);
            int top = CHUNK(v) == tc;
            if (bdir == 0) {
                int f = bi, k = imin(f, pending);
                bop(&b, RR, k); pending -= k; bop(&b, RA, f - k);
                if (top) { bop(&b, RB, pending); pending = 0; }
            } else {
                bop(&b, RRA, (sa - bi) % sa);
                if (top) { bop(&b, RB, pending); pending = 0; }
            }
            bop(&b, PB, 1);
            if (!top) pending++;
        }
    }
    if (g_pair) {  /* swap the remaining pairs (nearest first) so that A is cyclically sorted */
        for (;;) {
            int sa = b.s.a.len, bi = -1, bd = INT_MAX;
            for (int i = 0; i < sa; i++) {
                int x = dq_get(&b.s.a, i);
                if (g_pair[x] >= 0 && !b.swapped[x]) { int d = imin(i, sa - i); if (d < bd) { bd = d; bi = i; } }
            }
            if (bi < 0) break;
            int x = dq_get(&b.s.a, bi);
            if (bi <= sa - bi) { for (int k = 0; k < bi && !b.swapped[x]; k++) bop(&b, RA, 1); }
            else { for (int k = 0; k < sa - bi && !b.swapped[x]; k++) bop(&b, RRA, 1); }
            if (!b.swapped[x]) { b.swapped[x] = 1; sim_apply(&b.s, SA); ov_push(b.ops, SA); }
        }
    }
    if (bw > 0 && n <= 512) {
        int na = b.s.a.len, nb = b.s.b.len;
        int *Av = xmalloc(sizeof(int) * (na + 1)), *Bv = xmalloc(sizeof(int) * (nb + 1));
        for (int i = 0; i < na; i++) Av[i] = dq_get(&b.s.a, i);
        for (int i = 0; i < nb; i++) Bv[i] = dq_get(&b.s.b, i);
        OpVec rest; ov_init(&rest);
        beam_solve(n, Av, na, na ? Av[0] : 0, Bv, nb, bw, bc, &rest);
        for (int i = 0; i < rest.n; i++) bop(&b, rest.d[i], 1);
        ov_free(&rest); free(Av); free(Bv); sim_free(&b.s); free(crank); free(b.swapped);
        return;
    }
    int *idxA = xmalloc(sizeof(int) * n), *succ = xmalloc(sizeof(int) * n);
    char *inA = xmalloc(n);
    while (b.s.b.len) {
        int sa = b.s.a.len, sb = b.s.b.len;
        memset(inA, 0, n);
        for (int i = 0; i < sa; i++) { int x = dq_get(&b.s.a, i); idxA[x] = i; inA[x] = 1; }
        int nxt = -1;
        for (int v = n - 1; v >= 0; v--) { succ[v] = nxt; if (inA[v]) nxt = v; }
        int minA = nxt, bestc = INT_MAX, bj = -1, bfa = 0;
        for (int j = 0; j < sb; j++) {
            int v = dq_get(&b.s.b, j);
            int s2 = succ[v] >= 0 ? succ[v] : minA;
            int fa = sa ? idxA[s2] : 0;
            int c = comb_cost(fa, sa, j, sb, NULL);
            if (c < bestc) { bestc = c; bj = j; bfa = fa; }
        }
        int fa = bfa, fb = bj;
        int ra = fa, rra = fa ? sa - fa : 0, rb = fb, rrb = fb ? sb - fb : 0;
        int c0 = imax(ra, rb), c1 = imax(rra, rrb), c2 = ra + rrb, c3 = rra + rb;
        int mc = imin(imin(c0, c1), imin(c2, c3));
        if (mc == c0) { int k = imin(ra, rb); bop(&b, RR, k); bop(&b, RA, ra - k); bop(&b, RB, rb - k); }
        else if (mc == c1) { int k = imin(rra, rrb); bop(&b, RRR, k); bop(&b, RRA, rra - k); bop(&b, RRB, rrb - k); }
        else if (mc == c2) { bop(&b, RA, ra); bop(&b, RRB, rrb); }
        else { bop(&b, RRA, rra); bop(&b, RB, rb); }
        bop(&b, PA, 1);
    }
    int sa = b.s.a.len, i0 = 0;
    for (int i = 0; i < sa; i++) if (dq_get(&b.s.a, i) == 0) i0 = i;
    if (i0 <= sa - i0) bop(&b, RA, i0); else bop(&b, RRA, sa - i0);
    free(idxA); free(succ); free(inA); sim_free(&b.s); free(crank); free(b.swapped);
#undef CHUNK
}
