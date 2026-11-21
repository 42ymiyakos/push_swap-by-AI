/* Peephole optimizer: replace ops[i, i+K) by a shorter equivalent found by meet-in-the-middle BFS
 * over local stack windows (front/back parts; the middle block cannot move within K ops). */
#include "ps.h"

#define LCAP 64
typedef struct { int16_t f[LCAP], b[LCAP]; int nf, nb, mid; } LStack;
typedef struct { LStack A, B; } LState;

static int popf(LStack *s, int16_t *x) {
    if (s->nf) { *x = s->f[0]; memmove(s->f, s->f + 1, sizeof(int16_t) * (s->nf - 1)); s->nf--; return 1; }
    if (s->mid == 0 && s->nb) { *x = s->b[0]; memmove(s->b, s->b + 1, sizeof(int16_t) * (s->nb - 1)); s->nb--; return 1; }
    return 0;
}
static int popb(LStack *s, int16_t *x) {
    if (s->nb) { *x = s->b[--s->nb]; return 1; }
    if (s->mid == 0 && s->nf) { *x = s->f[--s->nf]; return 1; }
    return 0;
}
static int pushf(LStack *s, int16_t x) {
    if (s->nf >= LCAP) return 0;
    memmove(s->f + 1, s->f, sizeof(int16_t) * s->nf); s->f[0] = x; s->nf++; return 1;
}
static int pushb(LStack *s, int16_t x) { if (s->nb >= LCAP) return 0; s->b[s->nb++] = x; return 1; }
static int lsize(const LStack *s) { return s->nf + s->nb + s->mid; }
static int sw(LStack *s) {
    if (lsize(s) < 2) return 1;
    int16_t x, y;
    if (!popf(s, &x) || !popf(s, &y)) return 0;
    return pushf(s, x) && pushf(s, y);
}
static int rot(LStack *s) { if (lsize(s) < 2) return 1; int16_t x; if (!popf(s, &x)) return 0; return pushb(s, x); }
static int rrot(LStack *s) { if (lsize(s) < 2) return 1; int16_t x; if (!popb(s, &x)) return 0; return pushf(s, x); }
static int lapply(LState *st, int o) {
    int16_t x;
    switch (o) {
    case SA: return sw(&st->A);
    case SB: return sw(&st->B);
    case SS: return sw(&st->A) && sw(&st->B);
    case PA: if (!lsize(&st->B)) return 1; if (!popf(&st->B, &x)) return 0; return pushf(&st->A, x);
    case PB: if (!lsize(&st->A)) return 1; if (!popf(&st->A, &x)) return 0; return pushf(&st->B, x);
    case RA: return rot(&st->A);
    case RB: return rot(&st->B);
    case RR: return rot(&st->A) && rot(&st->B);
    case RRA: return rrot(&st->A);
    case RRB: return rrot(&st->B);
    case RRR: return rrot(&st->A) && rrot(&st->B);
    }
    return 0;
}
static const int INV[NOPS] = {SA, SB, SS, PB, PA, RRA, RRB, RRR, RA, RB, RR};

static uint64_t hs(const LState *s) {
    uint64_t h = 1469598103934665603ULL;
#define MIX(x) (h ^= (uint64_t)(x) + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2))
    const LStack *sts[2] = {&s->A, &s->B};
    for (int k = 0; k < 2; k++) {
        const LStack *st = sts[k];
        MIX(0xABCD); MIX((uint64_t)(int64_t)st->mid);
        if (st->mid == 0) { for (int i = 0; i < st->nf; i++) MIX(st->f[i]); for (int i = 0; i < st->nb; i++) MIX(st->b[i]); }
        else { for (int i = 0; i < st->nf; i++) MIX(st->f[i]); MIX(0x5555); for (int i = 0; i < st->nb; i++) MIX(st->b[i]); }
    }
#undef MIX
    return h;
}

static void mk(const Dq *d, int L, LStack *s) {
    int n = d->len;
    s->nf = s->nb = s->mid = 0;
    if (n <= 2 * L + 2) { for (int i = 0; i < n; i++) s->f[s->nf++] = (int16_t)dq_get(d, i); return; }
    for (int i = 0; i < L; i++) s->f[s->nf++] = (int16_t)dq_get(d, i);
    for (int i = n - L; i < n; i++) s->b[s->nb++] = (int16_t)dq_get(d, i);
    s->mid = n - 2 * L;
}

/* hash -> path table with insertion order */
typedef struct { uint64_t h; int idx; unsigned gen; } HSlot;
typedef struct {
    HSlot *tab; size_t mask; unsigned gen;
    uint64_t *ord; uint8_t (*path)[8]; uint8_t *plen; int n, cap;
} PathMap;
static void pm_init(PathMap *m, int cap) {
    size_t tc = 1; while (tc < (size_t)cap * 2) tc <<= 1;
    m->tab = xcalloc(tc, sizeof(HSlot)); m->mask = tc - 1; m->gen = 0;
    m->ord = xmalloc(sizeof(uint64_t) * cap); m->path = xmalloc(8 * (size_t)cap); m->plen = xmalloc(cap); m->cap = cap; m->n = 0;
}
static void pm_free(PathMap *m) { free(m->tab); free(m->ord); free(m->path); free(m->plen); }
static void pm_clear(PathMap *m) { m->gen++; m->n = 0; }
static int pm_find(const PathMap *m, uint64_t h) {
    size_t i = (size_t)(h * 0x9E3779B97F4A7C15ULL >> 32) & m->mask;
    while (m->tab[i].gen == m->gen) { if (m->tab[i].h == h) return m->tab[i].idx; i = (i + 1) & m->mask; }
    return -1;
}
static int pm_add(PathMap *m, uint64_t h, const uint8_t *p, int len) {
    if (m->n >= m->cap) return -1;
    size_t i = (size_t)(h * 0x9E3779B97F4A7C15ULL >> 32) & m->mask;
    while (m->tab[i].gen == m->gen) i = (i + 1) & m->mask;
    m->tab[i].gen = m->gen; m->tab[i].h = h; m->tab[i].idx = m->n;
    m->ord[m->n] = h; if (len) memcpy(m->path[m->n], p, len); m->plen[m->n] = (uint8_t)len;
    return m->n++;
}

typedef struct { LState s; uint8_t p[8]; int len; } BNode;

static void bfs(const LState *s0, int depth, PathMap *seen, BNode **cur, BNode **nxt, int lvcap) {
    pm_clear(seen);
    pm_add(seen, hs(s0), NULL, 0);
    int nc = 1, nn;
    (*cur)[0].s = *s0; (*cur)[0].len = 0;
    for (int d = 0; d < depth; d++) {
        nn = 0;
        for (int c = 0; c < nc; c++) for (int o = 0; o < NOPS; o++) {
            BNode *t = &(*nxt)[nn];
            t->s = (*cur)[c].s;
            if (!lapply(&t->s, o)) continue;
            uint64_t h = hs(&t->s);
            if (pm_find(seen, h) >= 0) continue;
            memcpy(t->p, (*cur)[c].p, (*cur)[c].len); t->p[(*cur)[c].len] = (uint8_t)o; t->len = (*cur)[c].len + 1;
            if (pm_add(seen, h, t->p, t->len) < 0) continue;
            if (nn + 1 < lvcap) nn++;
        }
        BNode *tmp = *cur; *cur = *nxt; *nxt = tmp; nc = nn;
    }
}

int peep_optimize(const int *init, int n, OpVec *ops, int K) {
    if (K > 10) K = 10;
    int saved = 0, d1 = K / 2, d2 = (K - 1) - d1;
    int cap = 1; for (int d = 0; d < d1 + 1; d++) cap *= NOPS; cap += 16;
    PathMap F, B; pm_init(&F, cap); pm_init(&B, cap);
    BNode *cur = xmalloc(sizeof(BNode) * cap), *nxt = xmalloc(sizeof(BNode) * cap);
    OpVec cand; ov_init(&cand);
    Sim s; sim_init(&s, init, n);
    for (int i = 0; i + K <= ops->n; sim_apply(&s, ops->d[i]), i++) {
        LState si, sj;
        mk(&s.a, K + 1, &si.A); mk(&s.b, K + 1, &si.B);
        sj = si;
        int ok = 1;
        for (int k = 0; k < K && ok; k++) ok = lapply(&sj, ops->d[i + k]);
        if (!ok) continue;
        bfs(&si, d1, &F, &cur, &nxt, cap);
        bfs(&sj, d2, &B, &cur, &nxt, cap);
        int best = K, bf = -1, bb = -1;
        for (int q = 0; q < B.n; q++) {
            int fi = pm_find(&F, B.ord[q]);
            if (fi < 0) continue;
            int len = F.plen[fi] + B.plen[q];
            if (len < best) { best = len; bf = fi; bb = q; }
        }
        if (best < K) {
            cand.n = 0;
            for (int k = 0; k < i; k++) ov_push(&cand, ops->d[k]);
            for (int k = 0; k < F.plen[bf]; k++) ov_push(&cand, F.path[bf][k]);
            for (int k = B.plen[bb] - 1; k >= 0; k--) ov_push(&cand, INV[B.path[bb][k]]);
            for (int k = i + K; k < ops->n; k++) ov_push(&cand, ops->d[k]);
            if (verify_ops(init, n, &cand)) { saved += ops->n - cand.n; ov_copy(ops, &cand); }
            if (i + K > ops->n) break;
        }
    }
    sim_free(&s); ov_free(&cand); free(cur); free(nxt); pm_free(&F); pm_free(&B);
    return saved;
}
