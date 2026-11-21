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

typedef struct { Sim s; OpVec *ops; } Builder;
static void bop(Builder *b, int o, int k) { for (int i = 0; i < k; i++) { sim_apply(&b->s, o); ov_push(b->ops, o); } }

void chunk_construct(const int *a0, int n, const char *inK, int m, int bw, int bc, OpVec *out) {
    Builder b; sim_init(&b.s, a0, n); b.ops = out; out->n = 0;
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
                int c = m - 1 - (int)((long long)v * m / n);
                if (c != tc && c != bcn) continue;
                int top = c == tc;
                int f = i, cf = top ? imax(f, pending) : f;
                int rv = (sa - i) % sa, cr = rv + (top ? pending : 0);
                int d = imin(cf, cr);
                if (d < bd) { bd = d; bi = i; bdir = cf <= cr ? 0 : 1; }
            }
            if (bi < 0) break;
            int v = dq_get(&b.s.a, bi);
            int top = (m - 1 - (int)((long long)v * m / n)) == tc;
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
    if (bw > 0 && n <= 512) {
        int na = b.s.a.len, nb = b.s.b.len;
        int *Av = xmalloc(sizeof(int) * (na + 1)), *Bv = xmalloc(sizeof(int) * (nb + 1));
        for (int i = 0; i < na; i++) Av[i] = dq_get(&b.s.a, i);
        for (int i = 0; i < nb; i++) Bv[i] = dq_get(&b.s.b, i);
        OpVec rest; ov_init(&rest);
        beam_solve(n, Av, na, na ? Av[0] : 0, Bv, nb, bw, bc, &rest);
        for (int i = 0; i < rest.n; i++) bop(&b, rest.d[i], 1);
        ov_free(&rest); free(Av); free(Bv); sim_free(&b.s);
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
    free(idxA); free(succ); free(inA); sim_free(&b.s);
}
