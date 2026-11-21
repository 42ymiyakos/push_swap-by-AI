/* Beam search for the B->A phase when A holds only settled (cyclically sorted) values. n <= 512. */
#include "ps.h"

#define BW 8
typedef struct { uint64_t w[BW]; } BS;

static inline void setb(BS *b, int i) { b->w[i >> 6] |= 1ULL << (i & 63); }
static inline void clrb(BS *b, int i) { b->w[i >> 6] &= ~(1ULL << (i & 63)); }
static inline int getb(const BS *b, int i) { return (int)(b->w[i >> 6] >> (i & 63) & 1); }
static int cnt(const BS *b, int l, int r) {
    if (l >= r) return 0;
    int c = 0, wl = l >> 6, wr = (r - 1) >> 6;
    if (wl == wr) return __builtin_popcountll(b->w[wl] >> (l & 63) & (~0ULL >> (63 - ((r - 1) & 63) + (l & 63))));
    c += __builtin_popcountll(b->w[wl] >> (l & 63));
    for (int w = wl + 1; w < wr; w++) c += __builtin_popcountll(b->w[w]);
    c += __builtin_popcountll(b->w[wr] & (~0ULL >> (63 - ((r - 1) & 63))));
    return c;
}
static int ccnt(const BS *b, int l, int r, int n) { return l <= r ? cnt(b, l, r) : cnt(b, l, n) + cnt(b, 0, r); }
static int next_set(const BS *b, int i, int n) {
    for (int k = 0; k < 2; k++) {
        int lo = k == 0 ? i : 0, hi = k == 0 ? n : i;
        for (int x = lo; x < hi;) {
            int w = x >> 6; uint64_t m = b->w[w] >> (x & 63);
            if (m) { int y = x + __builtin_ctzll(m); return y < hi ? y : -1; }
            x = (w + 1) << 6;
        }
    }
    return -1;
}
static int prev_set(const BS *b, int i, int n) {
    for (int k = 0; k < n; k++) { int x = ((i - k) % n + n) % n; if (getb(b, x)) return x; }
    return -1;
}

typedef struct { BS a, b; int ta, tb, cost, nA, nB, parent; int16_t v; uint8_t mode; int fa, fb; float ev; } Node;
typedef struct { float ev; int node; } Cand;
typedef struct { int c, pos, fa, fb, mode; } Ent;

static int cmp_int(const void *x, const void *y) { int a = *(const int *)x, b = *(const int *)y; return (a > b) - (a < b); }
static int cmp_ent(const void *x, const void *y) {
    const Ent *a = x, *b = y;
    if (a->c != b->c) return a->c < b->c ? -1 : 1;
    return (a->pos > b->pos) - (a->pos < b->pos);
}
static int cmp_cand(const void *x, const void *y) {
    const Cand *a = x, *b = y;
    if (a->ev != b->ev) return a->ev < b->ev ? -1 : 1;
    return (a->node > b->node) - (a->node < b->node);
}

static float heur(const Node *s, int n, int m, const int *Bvals, int *pts) {
    int L = s->nA;
    if (L <= 1) return 0.0f;
    BS bv; memset(&bv, 0, sizeof bv);
    for (int i = 0; i < m; i++) if (getb(&s->b, i)) setb(&bv, Bvals[i]);
    int np = 0;
    int first = next_set(&s->a, 0, n), prevv = -1, rk = 0, st = 0, e = 0;
    int lastv = prev_set(&s->a, n - 1, n);
    for (int x = first; x >= 0;) {
        if (x == s->ta) st = rk;
        int lo = prevv < 0 ? lastv : prevv;
        int mk = lo < x ? cnt(&bv, lo + 1, x) > 0 : (cnt(&bv, lo + 1, n) + cnt(&bv, 0, x)) > 0;
        if (mk) pts[np++] = rk;
        prevv = x; rk++;
        int nx = x + 1 < n ? next_set(&s->a, x + 1, n) : -1;
        if (nx <= x) break;
        x = nx;
    }
    if (np == 0) { int d = st; return (float)imin(d, L - d); }
    pts[np++] = st; pts[np++] = e;
    qsort(pts, np, sizeof(int), cmp_int);
    int k = 0;
    for (int i = 0; i < np; i++) if (k == 0 || pts[i] != pts[k - 1]) pts[k++] = pts[i];
    int gmax = -1, gi = 0;
    for (int i = 0; i < k; i++) {
        int g = (i + 1 < k ? pts[i + 1] : pts[0] + L) - pts[i];
        if (g > gmax) { gmax = g; gi = i; }
    }
    int l = pts[(gi + 1) % k], arc = L - gmax;
    int os = ((st - l) % L + L) % L, oe = ((e - l) % L + L) % L;
    int c1 = os + arc + (arc - oe), c2 = (arc - os) + arc + oe;
    return (float)imin(c1, c2);
}

static void put(OpVec *ops, int o, int k) { for (int i = 0; i < k; i++) ov_push(ops, o); }

void beam_solve(int n, const int *Avals, int na, int topA, const int *Bvals, int m, int width, int cand, OpVec *out) {
    int allCap = 1024, nall = 0;
    Node *all = xmalloc(sizeof(Node) * allCap);
    Node root; memset(&root, 0, sizeof root);
    for (int i = 0; i < na; i++) setb(&root.a, Avals[i]);
    for (int i = 0; i < m; i++) setb(&root.b, i);
    root.ta = topA; root.tb = 0; root.cost = 0; root.nA = na; root.nB = m; root.parent = -1;
    all[nall++] = root;
    int *layer = xmalloc(sizeof(int) * (width + 1)), nl = 1, *nxt = xmalloc(sizeof(int) * (width + 1));
    layer[0] = 0;
    int *pts = xmalloc(sizeof(int) * (n + 4));
    int *seen = xmalloc(sizeof(int) * (2 * m + 8));
    Ent *es = xmalloc(sizeof(Ent) * (2 * m + 8));
    int ccap = 1024, nc = 0;
    Cand *cands = xmalloc(sizeof(Cand) * ccap);
    uint64_t *hs = xmalloc(sizeof(uint64_t) * (width + 1));
    for (int depth = 0; depth < m; depth++) {
        nc = 0;
        for (int li = 0; li < nl; li++) {
            int id = layer[li];
            Node s = all[id];
            int sa = s.nA, sb = s.nB;
            int lim = imin(sb, 2 * cand + 40), ns = 0;
            for (int x = s.tb, k = 0; k < sb && k <= lim; k++) {
                x = next_set(&s.b, x, m); if (x < 0) break;
                seen[ns++] = x; x = (x + 1) % m;
            }
            for (int x = s.tb, k = 0; k < sb && k <= lim; k++) {
                x = prev_set(&s.b, (x - 1 + m) % m, m); if (x < 0) break;
                seen[ns++] = x;
            }
            qsort(seen, ns, sizeof(int), cmp_int);
            int u = 0;
            for (int i = 0; i < ns; i++) if (u == 0 || seen[i] != seen[u - 1]) seen[u++] = seen[i];
            int ne = 0;
            for (int q = 0; q < u; q++) {
                int pos = seen[q], v = Bvals[pos];
                int fb = ccnt(&s.b, s.tb, pos, m);
                int sv = next_set(&s.a, v, n);
                int fa = ccnt(&s.a, s.ta, sv, n);
                int md, c = comb_cost(fa, sa, fb, sb, &md);
                es[ne].c = c; es[ne].pos = pos; es[ne].fa = fa; es[ne].fb = fb; es[ne].mode = md; ne++;
            }
            qsort(es, ne, sizeof(Ent), cmp_ent);
            int take = 0;
            for (int k = 0; k < ne; k++) {
                if (take >= cand && es[k].c > es[take - 1].c) break;
                take++;
                Ent *e = &es[k];
                Node ch = s;
                int v = Bvals[e->pos];
                clrb(&ch.b, e->pos); ch.nB--;
                setb(&ch.a, v); ch.nA++; ch.ta = v;
                ch.tb = ch.nB ? next_set(&ch.b, (e->pos + 1) % m, m) : 0;
                ch.cost = s.cost + e->c + 1; ch.parent = id; ch.v = (int16_t)v; ch.mode = (uint8_t)e->mode; ch.fa = e->fa; ch.fb = e->fb;
                ch.ev = ch.cost + heur(&ch, n, m, Bvals, pts);
                if (nall == allCap) { allCap *= 2; all = xrealloc(all, sizeof(Node) * allCap); }
                all[nall++] = ch;
                if (nc == ccap) { ccap *= 2; cands = xrealloc(cands, sizeof(Cand) * ccap); }
                cands[nc].ev = ch.ev; cands[nc].node = nall - 1; nc++;
            }
        }
        qsort(cands, nc, sizeof(Cand), cmp_cand);
        int nn = 0, nh = 0;
        for (int c = 0; c < nc; c++) {
            const Node *s = &all[cands[c].node];
            uint64_t h = 1469598103934665603ULL;
            for (int w = 0; w < BW; w++) { h ^= s->b.w[w]; h *= 1099511628211ULL; }
            h ^= (uint64_t)s->ta * 1315423911u + s->tb; h *= 1099511628211ULL;
            int dup = 0;
            for (int q = 0; q < nh; q++) if (hs[q] == h) { dup = 1; break; }
            if (dup) continue;
            hs[nh++] = h; nxt[nn++] = cands[c].node;
            if (nn >= width) break;
        }
        int *t = layer; layer = nxt; nxt = t; nl = nn;
    }
    int bestId = -1, bestC = INT_MAX, bestF = 0;
    for (int li = 0; li < nl; li++) {
        const Node *s = &all[layer[li]];
        int f = cnt(&s->a, 0, s->ta), L = s->nA, fr = imin(f, (L - f) % L);
        if (s->cost + fr < bestC) { bestC = s->cost + fr; bestId = layer[li]; bestF = f; }
    }
    int plen = 0;
    for (int id = bestId; id > 0; id = all[id].parent) plen++;
    int *path = xmalloc(sizeof(int) * (plen + 1));
    { int k = plen; for (int id = bestId; id > 0; id = all[id].parent) path[--k] = id; }
    out->n = 0;
    for (int k = 0; k < plen; k++) {
        const Node *c = &all[path[k]], *p = &all[c->parent];
        int sa = p->nA, sb = p->nB;
        int ra = c->fa, rra = c->fa ? sa - c->fa : 0, rb = c->fb, rrb = c->fb ? sb - c->fb : 0;
        switch (c->mode) {
        case 0: put(out, RR, imin(ra, rb)); put(out, RA, ra - imin(ra, rb)); put(out, RB, rb - imin(ra, rb)); break;
        case 1: put(out, RRR, imin(rra, rrb)); put(out, RRA, rra - imin(rra, rrb)); put(out, RRB, rrb - imin(rra, rrb)); break;
        case 2: put(out, RA, ra); put(out, RRB, rrb); break;
        case 3: put(out, RRA, rra); put(out, RB, rb); break;
        }
        ov_push(out, PA);
    }
    int L = all[bestId].nA;
    if (bestF <= L - bestF) put(out, RRA, bestF); else put(out, RA, L - bestF);
    free(all); free(layer); free(nxt); free(pts); free(seen); free(es); free(cands); free(hs); free(path);
}
