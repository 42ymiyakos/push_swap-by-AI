/* Iterated Greedy with exact O(Q) DP reinsertion of one value (see ALGORITHM_NOTES.md).
 * Reinsertion of v chooses i (A2B before old query i), j (B2A before old query j, j >= i),
 * B slot b and A target slot t, or "keep v in A". Port of src/ig.hpp + IG in src/psig.cpp
 * (default configuration: all B slots via sweep, query-destroy, greedy acceptance). */
#include "ps.h"

typedef struct {
    int sa, sb;
    ll ancA, ancB;
    int lbA, lbB;
    ll src, dst;
    int srcIdx, dstRaw, t, c;
    int pr, prr;  /* pending A2A rotations before this query */
    int u;        /* moved value */
} QInfo;

typedef struct {
    QInfo *qi; int qcap;
    int *rkI, *P, *G, *Aarr, *bFbC2, *bFbS2, *bC1ii; int acap;
    int *rb; int rbcap;    /* nb x (Q+1) */
    int *rt; int rtcap;    /* nt x (Q+1) */
    int *tSUF, *tF2, *tFS, *tFC; int tcap;  /* nt x (Q+2) */
    int *tFinA, *tFinB; int tfcap;
    int *swR, *swG, *swPm, *swPmi, *slotBest, *slotI, *slotJ; int swcap;
    int Q, finSa, finTa, finOld, finLb, minIdx;
    ll finAnc;
    int nb, nt, W;  /* W = Q+1 row stride, Wt = Q+2 */
    int tmode, t0, topb, ubInit;
    int *tb; int tbcap;
    int *Dg, *Dr; int dcap;
    int *pbuf; size_t pbcap;  /* deferred range-add buffers for swG / swR */
    int allowA2A, a2aMax, irot;
} RI;

/* exact ((x % m) + m) % m, with a fast path for the common range (-m, 2m) */
static inline int md(int x, int m) {
    if (m <= 0) return 0;
    if (x < 0) { x += m; if (x < 0) x = ((x % m) + m) % m; }
    else if (x >= m) { x -= m; if (x >= m) x %= m; }
    return x;
}
/* comb_sync equals comb_cost when there are no pending A2A rotations */
static inline int combs(int fa, int sa, int fb, int sb, int pr, int prr) {
    return (pr | prr) ? comb_sync(fa, sa, fb, sb, pr, prr, NULL, NULL) : comb_cost(fa, sa, fb, sb, NULL);
}
static inline int mq(int x, int m) { x = x < 0 ? x + m : x; x = x >= m ? x - m : x; return m > 0 ? x : 0; }

static void grow(int **p, int *cap, int need) { if (*cap < need) { *cap = need * 2; *p = xrealloc(*p, sizeof(int) * *cap); } }

static void pass(RI *R, const Rep *r, ll EkI, const ll *Eb, int nb, const ll *Et, int nt, St *s) {
    int Q = r->nq, W = Q + 1;
    R->Q = Q; R->W = W; R->nb = nb; R->nt = nt;
    if (R->qcap < Q + 1) { R->qcap = 2 * (Q + 1); R->qi = xrealloc(R->qi, sizeof(QInfo) * R->qcap); }
    if (R->acap < Q + 2) {
        R->acap = 2 * (Q + 2);
        int **arrs[] = {&R->rkI, &R->P, &R->G, &R->Aarr, &R->bFbC2, &R->bFbS2, &R->bC1ii};
        for (int k = 0; k < 7; k++) *arrs[k] = xrealloc(*arrs[k], sizeof(int) * R->acap);
    }
    grow(&R->rb, &R->rbcap, nb * W + 1);
    grow(&R->rt, &R->rtcap, nt * W + 1);
    int ebSorted = 1;
    for (int k = 1; k < nb; k++) if (Eb[k - 1] > Eb[k]) ebSorted = 0;
    initial_state(r, s);
    ll ancA = s->na ? s->A[s->ta] : 0, ancB = 0;
    int ancT = 0;  /* 0: top = first element >= ancA; 1: top = last element < ancA */
    for (int i = 0; i <= Q; i++) {
        int sa = s->na, sb = s->nb;
        R->rkI[i] = lb_ll(s->A, sa, EkI);
        if (ebSorted) {
            int p = 0;
            for (int k = 0; k < nb; k++) { while (p < sb && s->B[p] < Eb[k]) p++; R->rb[k * W + i] = p; }
        } else for (int k = 0; k < nb; k++) R->rb[k * W + i] = lb_ll(s->B, sb, Eb[k]);
        for (int k = 0; k < nt; k++) R->rt[k * W + i] = lb_ll(s->A, sa, Et[k]);
        if (i == Q) {
            R->finSa = sa; R->finTa = s->ta; R->finAnc = ancA;
            R->finLb = lb_ll(s->A, sa, ancA) - ancT;
            R->finOld = final_rot(s, r, NULL);
            R->minIdx = r->minv < r->n ? lb_ll(s->A, sa, r->tk[r->minv] * 1024 + r->minv) : 0;
            break;
        }
        QInfo *f = &R->qi[i];
        f->sa = sa; f->sb = sb; f->ancA = ancA; f->ancB = ancB; f->pr = s->pr; f->prr = s->prr;
        f->lbA = lb_ll(s->A, sa, ancA) - ancT;
        f->lbB = lb_ll(s->B, sb, ancB);
        Query qu = r->q[i]; int u = qu.v; f->t = qu.t; f->u = u;
        if (qu.t == 2) {
            f->src = r->ik[u] * 1024 + u; f->dst = r->tk[u] * 1024 + u;
            f->srcIdx = lb_ll(s->A, sa, f->src); f->dstRaw = lb_ll(s->A, sa, f->dst);
            int wasTop = s->na && s->A[s->ta] == f->src;
            A2AInfo in;
            f->c = do_step_a2a(s, r, qu, &in);
            if (in.topWo == -2 || f->sa <= 1) { if (wasTop) { ancA = f->dst; ancT = 0; } }
            else if (in.fwd) { ancA = f->dst; ancT = 1; }
            else { ancA = f->dst; ancT = 0; }
            continue;
        }
        if (qu.t == 0) {
            f->src = r->ik[u] * 1024 + u; f->dst = r->bk[u] * 1024 + u;
            f->srcIdx = lb_ll(s->A, sa, f->src); f->dstRaw = lb_ll(s->B, sb, f->dst);
        } else {
            f->src = r->bk[u] * 1024 + u; f->dst = r->tk[u] * 1024 + u;
            f->srcIdx = lb_ll(s->B, sb, f->src); f->dstRaw = lb_ll(s->A, sa, f->dst);
        }
        f->c = do_step(s, r, qu, NULL);
        if (qu.t == 0) { ancA = f->src; ancB = f->dst; } else { ancA = f->dst; ancB = f->src; }
        ancT = 0;
    }
}

/* cost of old query i with optional extra element in A (xa) and/or B (xb); top overrides (-1 = anchor) */
static inline int aug(const RI *R, int i, ll xa, ll xb, int topAo, int topBo) {
    const QInfo *f = &R->qi[i];
    int sa = f->sa + (xa >= 0), sb = f->sb + (xb >= 0);
    int ta = topAo >= 0 ? topAo : md(f->lbA + (xa >= 0 && xa < f->ancA), sa);
    if (f->t == 2) return a2a_cost(sa, ta, f->srcIdx + (xa >= 0 && xa < f->src), f->dstRaw + (xa >= 0 && xa < f->dst)).cost;
    int tb = topBo >= 0 ? topBo : md(f->lbB + (xb >= 0 && xb < f->ancB), sb);
    int fa, fb;
    if (f->t == 0) {
        int si = f->srcIdx + (xa >= 0 && xa < f->src);
        fa = md(si - ta, sa);
        int sl = sb ? md(f->dstRaw + (xb >= 0 && xb < f->dst), sb) : 0;
        fb = sb ? md(sl - tb, sb) : 0;
    } else {
        int si = f->srcIdx + (xb >= 0 && xb < f->src);
        fb = md(si - tb, sb);
        int sl = sa ? md(f->dstRaw + (xa >= 0 && xa < f->dst), sa) : 0;
        fa = sa ? md(sl - ta, sa) : 0;
    }
    return 1 + combs(fa, sa, fb, sb, f->pr, f->prr);
}
static inline int augA(const RI *R, int i, ll xa, int ta) {
    const QInfo *f = &R->qi[i];
    int sa = f->sa + 1;
    if (f->t == 0) return md(f->srcIdx + (xa < f->src) - ta, sa);
    return md(md(f->dstRaw + (xa < f->dst), sa) - ta, sa);
}
static inline int augB(const RI *R, int i, int tb) {
    const QInfo *f = &R->qi[i];
    int sb = f->sb;
    if (!sb) return 0;
    if (f->t == 0) return md(md(f->dstRaw, sb) - tb, sb);
    return md(f->srcIdx - tb, sb);
}

typedef struct { const RI *R; const Rep *r; int v; } FinCtx;
static int fin(const RI *R, const Rep *r, int v, ll x, int rx, int topo) {
    int sa = R->finSa + 1;
    int ta = topo >= 0 ? topo : md(R->finLb + (x < R->finAnc), sa);
    int tgt;
    if (v < r->minv) tgt = rx;
    else { ll em = r->tk[r->minv] * 1024 + r->minv; tgt = R->minIdx + (x < em); }
    int f = md(tgt - ta, sa);
    return imin(f, sa - f) - R->finOld;
}

/* per-slot arrays live in R->pbuf (5 arrays of W ints per slot) so that a slot is prepared only once */
#define PB(R, bi, k) ((R)->pbuf + ((size_t)(bi) * 5 + (k)) * (R)->W)
static void prep_b(RI *R, int bi, ll EkI, ll eb, int vbf) {
    int Q = R->Q, W = R->W;
    const int *Rr = R->rb + bi * W;
    int *G = PB(R, bi, 0), *Aarr = PB(R, bi, 1), *bFbC2 = PB(R, bi, 2), *bFbS2 = PB(R, bi, 3), *bC1ii = PB(R, bi, 4);
    G[0] = 0;
    for (int q = 0; q < Q; q++) G[q + 1] = G[q] + aug(R, q, -1, eb, -1, -1) - R->qi[q].c;
    for (int i = 0; i <= Q; i++) {
        int sa = i < Q ? R->qi[i].sa : R->finSa, sb = i < Q ? R->qi[i].sb : 0;
        int sa1 = sa + 1, ta1;
        if (i == 0 && vbf) ta1 = R->rkI[0];
        else ta1 = i < Q ? md(R->qi[i].lbA + (EkI < R->qi[i].ancA), sa1) : md(R->finLb + (EkI < R->finAnc), sa1);
        int fa = md(R->rkI[i] - ta1, sa1);
        int fb = (i < Q && sb) ? md(md(Rr[i], sb) - md(R->qi[i].lbB, sb), sb) : 0;
        int c1 = 1 + combs(fa, sa1, fb, sb, i < Q ? R->qi[i].pr : 0, i < Q ? R->qi[i].prr : 0);
        bC1ii[i] = c1;
        if (i < Q) {
            const QInfo *f = &R->qi[i];
            int s1 = aug(R, i, -1, eb, f->sa ? md(R->rkI[i], f->sa) : 0, Rr[i]);
            Aarr[i] = R->P[i] + c1 + s1 - f->c - G[i + 1];
            int sb1 = f->sb + 1;
            int tb1 = md(f->lbB + (eb < f->ancB), sb1);
            bFbC2[i] = md(Rr[i] - tb1, sb1);
            bFbS2[i] = augB(R, i, f->sb ? md(Rr[i], f->sb) : 0);
        } else { bFbC2[i] = 0; bFbS2[i] = 0; }
    }
}

static void dp_bt(const RI *R, int bi, int ti, Best *bst) {
    int Q = R->Q, Wt = Q + 2, W = R->W;
    const int *G = PB(R, bi, 0), *Aarr = PB(R, bi, 1), *bFbC2 = PB(R, bi, 2), *bFbS2 = PB(R, bi, 3), *bC1ii = PB(R, bi, 4);
    const int *SU = R->tSUF + ti * Wt, *F2 = R->tF2 + ti * W, *FS = R->tFS + ti * W, *FC = R->tFC + ti * W;
    int finA = R->tFinA[ti], finB = R->tFinB[ti];
    int pm = INT_MAX, pmi = -1;
    for (int j = 0; j <= Q; j++) {
        if (j >= 1 && Aarr[j - 1] < pm) { pm = Aarr[j - 1]; pmi = j - 1; }
        int tail, c2;
        if (j < Q) {
            const QInfo *f = &R->qi[j];
            c2 = 1 + combs(F2[j], f->sa, bFbC2[j], f->sb + 1, f->pr, f->prr);
            int s2 = f->t == 2 ? FS[j] : 1 + comb_cost(FS[j], f->sa + 1, bFbS2[j], f->sb, NULL);
            tail = s2 - f->c + SU[j + 1] + finA;
        } else {
            c2 = 1 + comb_cost(F2[Q], R->finSa, 0, 1, NULL);
            tail = finB;
        }
        int Bj = G[j] + c2 + tail;
        if (pmi >= 0 && pm + Bj < bst->delta) { bst->delta = pm + Bj; bst->i = pmi; bst->j = j; bst->bi = bi; bst->ti = ti; bst->keep = 0; bst->a2a = 0; }
        int sa = j < Q ? R->qi[j].sa : R->finSa, sb = j < Q ? R->qi[j].sb : 0;
        int cc = 1 + comb_cost(FC[j], sa, 0, sb + 1, NULL);
        int d = R->P[j] + bC1ii[j] + cc + tail;
        if (d < bst->delta) { bst->delta = d; bst->i = j; bst->j = j; bst->bi = bi; bst->ti = ti; bst->keep = 0; bst->a2a = 0; }
    }
}

typedef struct { int d, k; } Pair;
static int cmp_pair(const void *x, const void *y) {
    const Pair *a = x, *b = y;
    if (a->d != b->d) return a->d < b->d ? -1 : 1;
    return (a->k > b->k) - (a->k < b->k);
}

/* everything that does not depend on the B slots: query pass, prefix sums, keep/A2A options, per-target arrays */
static Best solve_prep(RI *R, const Rep *r, int v, const ll *Eb, int nb, const ll *Et, int nt, int keepOK, int vbf, St *s) {
    ll EkI = r->ik[v] * 1024 + v;
    pass(R, r, EkI, Eb, nb, Et, nt, s);
    int Q = R->Q, W = R->W, Wt = Q + 2;
    Best best = {INT_MAX, -1, -1, -1, -1, 0, 0};
    R->P[0] = 0;
    for (int i = 0; i < Q; i++) {
        int topo = (i == 0 && vbf) ? R->rkI[0] : -1;
        R->P[i + 1] = R->P[i] + aug(R, i, EkI, -1, topo, -1) - R->qi[i].c;
    }
    if (keepOK) {
        int topo = (Q == 0 && vbf) ? R->rkI[0] : -1;
        int d = R->P[Q] + fin(R, r, v, EkI, R->rkI[Q], topo);
        if (d < best.delta) { best.delta = d; best.keep = 1; }
    }
    grow(&R->tSUF, &R->tcap, nt * Wt + 1);
    if (R->tcap < nt * Wt + 1) R->tcap = nt * Wt + 1;
    R->tF2 = xrealloc(R->tF2, sizeof(int) * (nt * W + 1));
    R->tFS = xrealloc(R->tFS, sizeof(int) * (nt * W + 1));
    R->tFC = xrealloc(R->tFC, sizeof(int) * (nt * W + 1));
    grow(&R->tFinA, &R->tfcap, nt + 1);
    R->tFinB = xrealloc(R->tFinB, sizeof(int) * R->tfcap);
    for (int ti = 0; ti < nt; ti++) {
        ll et = Et[ti];
        const int *T = R->rt + ti * W;
        int *SU = R->tSUF + ti * Wt;
        SU[Q] = 0; SU[Q + 1] = 0;
        for (int q = Q - 1; q >= 0; q--) SU[q] = SU[q + 1] + aug(R, q, et, -1, -1, -1) - R->qi[q].c;
        R->tFinA[ti] = fin(R, r, v, et, T[Q], -1);
        R->tFinB[ti] = fin(R, r, v, et, T[Q], T[Q]);
        int *F2 = R->tF2 + ti * W, *FS = R->tFS + ti * W, *FC = R->tFC + ti * W;
        for (int j = 0; j <= Q; j++) {
            int sa = j < Q ? R->qi[j].sa : R->finSa;
            int ta = j < Q ? md(R->qi[j].lbA, sa) : md(R->finLb, R->finSa);
            F2[j] = sa ? md(md(T[j], sa) - ta, sa) : 0;
            FS[j] = j < Q ? (R->qi[j].t == 2 ? aug(R, j, et, -1, T[j], -1) : augA(R, j, et, T[j])) : 0;
            int tA = sa ? md(R->rkI[j], sa) : 0;
            FC[j] = sa ? md(md(T[j], sa) - tA, sa) : 0;
        }
    }
    if (R->allowA2A) {
        for (int ti = 0; ti < nt; ti++) {
            ll et = Et[ti];
            const int *T = R->rt + ti * W, *SU = R->tSUF + ti * Wt;
            for (int i = 0; i <= Q; i++) {
                int sa = i < Q ? R->qi[i].sa : R->finSa, sa1 = sa + 1;
                if (sa < 1) continue;
                int ta1 = (i == 0 && vbf) ? R->rkI[0] : (i < Q ? md(R->qi[i].lbA + (EkI < R->qi[i].ancA), sa1) : md(R->finLb + (EkI < R->finAnc), sa1));
                int src = R->rkI[i], dr = T[i] + (EkI < et);
                A2AInfo in = a2a_cost(sa1, ta1, src, dr);
                if (in.topWo == -2) continue;
                if (in.steps > R->a2aMax) continue;
                int nt2 = in.fwd ? in.topWo + (T[i] <= in.topWo) : T[i];
                int d;
                if (i < Q) {
                    int s3 = aug(R, i, et, -1, nt2, -1);
                    d = R->P[i] + in.cost + s3 - R->qi[i].c + SU[i + 1] + R->tFinA[ti];
                } else d = R->P[Q] + in.cost + fin(R, r, v, et, T[Q], nt2);
                if (d < best.delta) { best.delta = d; best.i = i; best.j = i; best.bi = -1; best.ti = ti; best.keep = 0; best.a2a = 1; }
            }
        }
    }
    return best;
}

/* best over the given B slots (ranks must be in R->rb), starting from a base candidate */
static Best solve_slots(RI *R, const ll *Eb, int nb, int nt, int vbf, ll EkI, Best best) {
    {   size_t need = (size_t)(nb + 1) * 5 * R->W;
        if (R->pbcap < need) { R->pbcap = 2 * need; R->pbuf = xrealloc(R->pbuf, sizeof(int) * R->pbcap); } }
    if (R->tmode == 0 || nt == 1) {
        for (int bi = 0; bi < nb; bi++) { prep_b(R, bi, EkI, Eb[bi], vbf); for (int ti = 0; ti < nt; ti++) dp_bt(R, bi, ti, &best); }
    } else {
        int tp = imin(imax(R->t0, 0), nt - 1);
        Pair *cand = xmalloc(sizeof(Pair) * (nb + 1));
        for (int bi = 0; bi < nb; bi++) {
            prep_b(R, bi, EkI, Eb[bi], vbf);
            Best bb = {INT_MAX, -1, -1, -1, -1, 0, 0};
            dp_bt(R, bi, tp, &bb);
            cand[bi].d = bb.delta; cand[bi].k = bi;
            if (bb.delta < best.delta) best = bb;
        }
        int k = imin(nb, R->topb);
        qsort(cand, nb, sizeof(Pair), cmp_pair);
        for (int c = 0; c < k; c++)  /* slot arrays are still in pbuf: no recomputation */
            for (int ti = 0; ti < nt; ti++) dp_bt(R, cand[c].k, ti, &best);
        free(cand);
    }
    return best;
}

static Best solve(RI *R, const Rep *r, int v, const ll *Eb, int nb, const ll *Et, int nt, int keepOK, int vbf, St *s) {
    Best base = solve_prep(R, r, v, Eb, nb, Et, nt, keepOK, vbf, s);
    return solve_slots(R, Eb, nb, nt, vbf, r->ik[v] * 1024 + v, base);
}

/* ranks of bl-slots (slot q lies after bl[q]) over time, from the query list: O(Q) per slot */
static void slot_ranks(RI *R, const int *slots, int k, const int *posOf) {
    int Q = R->Q, W = R->W;
    grow(&R->rb, &R->rbcap, k * W + 1);
    for (int z = 0; z < k; z++) {
        int q = slots[z], rk = 0;
        int *row = R->rb + z * W;
        for (int x = 0; x < Q; x++) {
            row[x] = rk;
            const QInfo *f = &R->qi[x];
            if (f->t == 0) { if (posOf[f->u] <= q) rk++; }
            else if (f->t == 1) { if (posOf[f->u] <= q) rk--; }
        }
        row[Q] = rk;
    }
}


/* all-slot sweep with provisional target ti0; slot s lies between bl[s] and bl[s+1]. */
static void sweep(RI *R, const Rep *r, int v, int m, const int *posOf, int ti0, int vbf) {
    ll EkI = r->ik[v] * 1024 + v;
    int Q = R->Q, W = R->W, Wt = Q + 2;
    const int *SU = R->tSUF + ti0 * Wt, *F2 = R->tF2 + ti0 * W, *FS = R->tFS + ti0 * W, *FC = R->tFC + ti0 * W;
    int finA = R->tFinA[ti0], finB = R->tFinB[ti0];
    if (R->swcap < m) {
        R->swcap = 2 * m;
        int **arrs[] = {&R->swR, &R->swG, &R->swPm, &R->swPmi, &R->slotBest, &R->slotI, &R->slotJ};
        for (int k = 0; k < 7; k++) *arrs[k] = xrealloc(*arrs[k], sizeof(int) * R->swcap);
    }
    int *swR = R->swR, *swG = R->swG, *swPm = R->swPm, *swPmi = R->swPmi, *slotBest = R->slotBest, *slotI = R->slotI, *slotJ = R->slotJ;
    for (int s = 0; s < m; s++) { swR[s] = 0; swG[s] = 0; swPm[s] = INT_MAX; swPmi[s] = -1; slotBest[s] = INT_MAX; slotI[s] = -1; slotJ[s] = -1; }
    int bestSoFar = R->ubInit;
    if (R->dcap < m + 1) { R->dcap = 2 * (m + 1); R->Dg = xrealloc(R->Dg, sizeof(int) * R->dcap); R->Dr = xrealloc(R->Dr, sizeof(int) * R->dcap); }
    int *Dg = R->Dg, *Dr = R->Dr, pendG = 0, pendR = 0;
    for (int s = 0; s <= m; s++) { Dg[s] = 0; Dr[s] = 0; }
    /* lbMin <= min over slots with a prefix of (swPm + swG); NOPREF when no slot has a prefix yet */
    const int NOPREF = INT_MAX / 4;
    int lbMin = NOPREF;
#define MATERIALIZE() do { \
        if (pendG) { int acc = 0; for (int s_ = 0; s_ < m; s_++) { acc += Dg[s_]; Dg[s_] = 0; swG[s_] += acc; } Dg[m] = 0; pendG = 0; } \
        if (pendR) { int acc = 0; for (int s_ = 0; s_ < m; s_++) { acc += Dr[s_]; Dr[s_] = 0; swR[s_] += acc; } Dr[m] = 0; pendR = 0; } \
    } while (0)
    for (int x = 0; x <= Q; x++) {
        int sa = x < Q ? R->qi[x].sa : R->finSa, sb = x < Q ? R->qi[x].sb : 0;
        int sa1 = sa + 1;
        int ta1 = (x == 0 && vbf) ? R->rkI[0] : (x < Q ? md(R->qi[x].lbA + (EkI < R->qi[x].ancA), sa1) : md(R->finLb + (EkI < R->finAnc), sa1));
        int fa1 = md(R->rkI[x] - ta1, sa1);
        int tbOld = x < Q ? md(R->qi[x].lbB, sb) : 0;
        int tA = sa ? md(R->rkI[x], sa) : 0;
        int cc = 1 + comb_cost(FC[x], sa, 0, sb + 1, NULL);
        /* A2B(v) at a time where A's top is far from v is not considered (irot >= 0) */
        int farA = R->irot >= 0 && imin(fa1, sa1 - fa1) > R->irot;
        if (x == Q) {
            MATERIALIZE();
            int c2 = 1 + comb_cost(F2[Q], R->finSa, 0, 1, NULL);
            for (int s = 0; s < m; s++) {
                int c1 = 1 + comb_cost(fa1, sa1, 0, 0, NULL);
                int d = R->P[Q] + c1 + cc + finB;
                if (!farA && d < slotBest[s]) { slotBest[s] = d; slotI[s] = Q; slotJ[s] = Q; }
                if (swPmi[s] >= 0) {
                    int dd = swPm[s] + swG[s] + c2 + finB;
                    if (dd < slotBest[s]) { slotBest[s] = dd; slotI[s] = swPmi[s]; slotJ[s] = Q; }
                }
            }
            break;
        }
        const QInfo *f = &R->qi[x];
        int sb1 = sb + 1;
        int pA = sb ? posOf[f->ancB & 1023] : -1;
        int u = (int)((f->t == 0 ? f->dst : f->src) & 1023);
        int pU = posOf[u];
        int faOld = f->t == 0 ? md(f->srcIdx - md(f->lbA, sa), sa) : (sa ? md(md(f->dstRaw, sa) - md(f->lbA, sa), sa) : 0);
        int faS1 = f->t == 0 ? md(f->srcIdx - tA, sa) : (sa ? md(md(f->dstRaw, sa) - tA, sa) : 0);
        int ev[2][2];
        for (int a = 0; a < 2; a++) for (int b2 = 0; b2 < 2; b2++) {
            int tb_ = sb ? md(f->lbB + a, sb1) : 0, fb;
            if (f->t == 0) { int sl = md(f->dstRaw + b2, sb1); fb = md(sl - tb_, sb1); }
            else { int si = f->srcIdx + b2; fb = md(si - tb_, sb1); }
            ev[a][b2] = f->t == 2 ? 0 : 1 + combs(faOld, sa, fb, sb1, f->pr, f->prr) - f->c;
        }
        int isAA = f->t == 2;
        int s1AA = isAA ? a2a_cost(sa, tA, f->srcIdx, f->dstRaw).cost : 0;
        int sync = f->pr | f->prr;
        int base = SU[x + 1] + finA - f->c;
        int lbS2 = isAA ? FS[x] : 1 + imin(FS[x], sa + 1 - FS[x]);
        int lbC2 = 1 + (sa ? imin(F2[x], sa - F2[x]) : 0) - (f->pr + f->prr);
        if (lbC2 < 1) lbC2 = 1;
        int lbTail = lbS2 + base;
        {   /* rank update of this query (deferred) */
            int dr = f->t == 0 ? 1 : (f->t == 1 ? -1 : 0);
            if (dr && pU >= 0 && pU < m) { Dr[pU] += dr; Dr[m] -= dr; pendR = 1; }
        }
        if (farA && (lbMin >= NOPREF || lbMin + lbC2 + lbTail >= bestSoFar)) {
            /* provably no slot passes the B2A test at this x and A2B(v) is not allowed here:
             * only the middle-effect prefix sums change; defer them as range adds */
            int lo = imin(pA, pU), hi = imax(pA, pU);
            lo = lo < 0 ? 0 : (lo > m ? m : lo);
            hi = hi < 0 ? 0 : (hi > m ? m : hi);
            int eLo = ev[1][1], eMid = pA > pU ? ev[1][0] : ev[0][1], eHi = ev[0][0];
            int emin = INT_MAX;
            if (lo > 0) { Dg[0] += eLo; Dg[lo] -= eLo; emin = imin(emin, eLo); }
            if (hi > lo) { Dg[lo] += eMid; Dg[hi] -= eMid; emin = imin(emin, eMid); }
            if (m > hi) { Dg[hi] += eHi; Dg[m] -= eHi; emin = imin(emin, eHi); }
            pendG = 1;
            if (lbMin < NOPREF && emin != INT_MAX) lbMin += emin;
            continue;
        }
        /* the rank update of this x must not be visible yet: undo its deferral, materialize, redo after */
        {
            int dr = f->t == 0 ? 1 : (f->t == 1 ? -1 : 0);
            if (dr && pU >= 0 && pU < m) { Dr[pU] -= dr; Dr[m] += dr; }
        }
        MATERIALIZE();
        /* c1 (needed for every slot) depends only on the slot's rank R in B: tabulate it (O(sb) <= O(m)) */
        int nR = sb + 1;
        if (R->tbcap < nR + 1) { R->tbcap = 2 * (nR + 1); R->tb = xrealloc(R->tb, sizeof(int) * R->tbcap); }
        int *Tc1 = R->tb;
        for (int Rv = 0; Rv < nR; Rv++) {
            int Rm = sb ? (Rv < sb ? Rv : Rv - sb) : 0;
            int fb1 = sb ? mq(Rm - tbOld, sb) : 0;
            Tc1[Rv] = 1 + (sync ? combs(fa1, sa1, fb1, sb, f->pr, f->prr) : comb_cost(fa1, sa1, fb1, sb, NULL));
        }
        if (farA) {
            for (int s = 0; s < m; s++) {
                int ltA = s < pA, ltU = s < pU;
                if (swPmi[s] >= 0 && swPm[s] + swG[s] + lbC2 + lbTail < bestSoFar) {
                    int Rv = swR[s];
                    int Rm = sb ? (Rv < sb ? Rv : Rv - sb) : 0;
                    int tb1 = sb ? mq(f->lbB + ltA, sb1) : 0;
                    int fbS2 = 0;
                    if (sb) fbS2 = f->t == 0 ? mq(mq(f->dstRaw, sb) - Rm, sb) : mq(f->srcIdx - Rm, sb);
                    int s2 = isAA ? FS[x] : 1 + comb_cost(FS[x], sa + 1, fbS2, sb, NULL);
                    int tail = s2 + base;
                    int fbC2 = mq(Rv - tb1, sb1);
                    int c2 = 1 + (sync ? combs(F2[x], sa, fbC2, sb1, f->pr, f->prr) : comb_cost(F2[x], sa, fbC2, sb1, NULL));
                    int dd = swPm[s] + swG[s] + c2 + tail;
                    if (dd < slotBest[s]) { slotBest[s] = dd; slotI[s] = swPmi[s]; slotJ[s] = x; }
                    if (slotBest[s] < bestSoFar) bestSoFar = slotBest[s];
                }
                swG[s] += ev[ltA][ltU];
            }
        } else
        for (int s = 0; s < m; s++) {
            int Rv = swR[s];
            int ltA = s < pA, ltU = s < pU;
            int Rm = sb ? (Rv < sb ? Rv : Rv - sb) : 0;
            int c1 = Tc1[Rv];
            int lbJ = imin(swPmi[s] >= 0 ? swPm[s] + swG[s] + lbC2 : INT_MAX / 2, R->P[x] + c1 + cc) + lbTail;
            if (lbJ < bestSoFar) {
                int tb1 = sb ? mq(f->lbB + ltA, sb1) : 0;
                int fbS2 = 0;
                if (sb) fbS2 = f->t == 0 ? mq(mq(f->dstRaw, sb) - Rm, sb) : mq(f->srcIdx - Rm, sb);
                int s2 = isAA ? FS[x] : 1 + comb_cost(FS[x], sa + 1, fbS2, sb, NULL);
                int tail = s2 + base;
                int d = R->P[x] + c1 + cc + tail;
                if (d < slotBest[s]) { slotBest[s] = d; slotI[s] = x; slotJ[s] = x; }
                if (swPmi[s] >= 0) {
                    int fbC2 = mq(Rv - tb1, sb1);
                    int c2 = 1 + (sync ? combs(F2[x], sa, fbC2, sb1, f->pr, f->prr) : comb_cost(F2[x], sa, fbC2, sb1, NULL));
                    int dd = swPm[s] + swG[s] + c2 + tail;
                    if (dd < slotBest[s]) { slotBest[s] = dd; slotI[s] = swPmi[s]; slotJ[s] = x; }
                }
                if (slotBest[s] < bestSoFar) bestSoFar = slotBest[s];
            }
            int e = ev[ltA][ltU];
            int pre = R->P[x] + c1 - f->c - (swG[s] + e);
            if (pre + 1 < swPm[s]) {
                int s1;
                if (isAA) s1 = s1AA;
                else {
                    int fbs1;
                    if (f->t == 0) { int sl = mq(f->dstRaw + ltU, sb1); fbs1 = mq(sl - Rv, sb1); }
                    else { int si = f->srcIdx + ltU; fbs1 = mq(si - Rv, sb1); }
                    s1 = 1 + comb_cost(faS1, sa, fbs1, sb1, NULL);
                }
                int Ax = pre + s1;
                if (Ax < swPm[s]) { swPm[s] = Ax; swPmi[s] = x; }
            }
            swG[s] += e;
        }
        {
            int dr = f->t == 0 ? 1 : (f->t == 1 ? -1 : 0);
            if (dr && pU >= 0 && pU < m) { Dr[pU] += dr; Dr[m] -= dr; pendR = 1; }
        }
        lbMin = NOPREF;
        for (int s = 0; s < m; s++) if (swPmi[s] >= 0) { int t = swPm[s] + swG[s]; if (t < lbMin) lbMin = t; }
    }
}

/* ---------------- IG driver ---------------- */
typedef struct {
    const Params *P; Rng *rng; RI R; St st;
    int *bl, *pos, *vcost, *D, *dc, *ord; char *mark;
    ll *Eb, *Et, *Eub, *Eb2;
    struct AE { ll k; int v, kind; } *al;
    int n;
} IG;

static int between(ll a, ll x, ll b) { return a < b ? (a < x && x < b) : (x > a || x < b); }

static const Rep *g_r;  /* comparator context */
static int cmp_bl(const void *x, const void *y) {
    int a = *(const int *)x, b = *(const int *)y;
    ll ka = g_r->bk[a], kb = g_r->bk[b];
    if (ka != kb) return ka < kb ? -1 : 1;
    return (a > b) - (a < b);
}
static int cmp_ae3(const void *x, const void *y) {
    const struct AE *a = x, *b = y;
    if (a->k != b->k) return a->k < b->k ? -1 : 1;
    if (a->v != b->v) return a->v < b->v ? -1 : 1;
    return (a->kind > b->kind) - (a->kind < b->kind);
}
static int cmp_ae2(const void *x, const void *y) {
    const struct AE *a = x, *b = y;
    if (a->k != b->k) return a->k < b->k ? -1 : 1;
    return (a->v > b->v) - (a->v < b->v);
}

static int build_bl(IG *g, const Rep *r) {
    int m = 0;
    for (int u = 0; u < r->n; u++) if (r->present[u] && !r->inK[u] && !r->isA2A[u]) g->bl[m++] = u;
    g_r = r; qsort(g->bl, m, sizeof(int), cmp_bl);
    return m;
}

static void respace(IG *g, Rep *r) {
    int m = build_bl(g, r);
    for (int i = 0; i < m; i++) r->bk[g->bl[i]] = KEY_BASE + (ll)i * KEY_GAP;
    int L = 0;
    for (int v = 0; v < r->n; v++) {
        int both = r->present[v] && r->inK[v];
        g->al[L].k = r->ik[v]; g->al[L].v = v; g->al[L].kind = both ? 2 : 0; L++;
        if (r->present[v] && !r->inK[v]) { g->al[L].k = r->tk[v]; g->al[L].v = v; g->al[L].kind = 1; L++; }
    }
    qsort(g->al, L, sizeof(struct AE), cmp_ae3);
    for (int i = 0; i < L; i++) {
        ll k = KEY_BASE + (ll)i * KEY_GAP;
        if (g->al[i].kind == 0) r->ik[g->al[i].v] = k;
        else if (g->al[i].kind == 1) r->tk[g->al[i].v] = k;
        else r->ik[g->al[i].v] = r->tk[g->al[i].v] = k;
    }
}

static void remove_value(Rep *r, int v) {
    if (!r->inK[v]) rep_remove_queries(r, v);
    r->present[v] = 0;
    rep_fix_minv(r);
}

static int reinsert(IG *g, Rep *r, int v) {
    RI *R = &g->R;
    for (int attempt = 0; attempt < 2; attempt++) {
        int m = build_bl(g, r);
        for (int u = 0; u < r->n; u++) g->pos[u] = -1;
        for (int k = 0; k < m; k++) g->pos[g->bl[k]] = k;
        int nEb = 0, bad = 0;
        int ns = m > 0 ? m : 1;
        for (int s = 0; s < ns; s++) {
            ll key;
            if (m == 0) key = KEY_BASE;
            else if (s < m - 1) { ll a = r->bk[g->bl[s]], b = r->bk[g->bl[s + 1]]; if (b - a < 2) { bad = 1; break; } key = a + (b - a) / 2; }
            else key = r->bk[g->bl[m - 1]] + KEY_GAP;
            g->Eb[nEb++] = key * 1024 + v;
        }
        int L = 0;
        for (int u = 0; u < r->n; u++) {
            g->al[L].k = r->ik[u]; g->al[L].v = u; g->al[L].kind = 0; L++;
            if (r->present[u] && !r->inK[u]) { g->al[L].k = r->tk[u]; g->al[L].v = u; g->al[L].kind = 0; L++; }
        }
        qsort(g->al, L, sizeof(struct AE), cmp_ae2);
        int pv = -1, sv = -1;
        for (int d = 1; d < r->n; d++) { int w = (v - d + r->n) % r->n; if (r->present[w] && !r->ph[w]) { pv = w; break; } }
        for (int d = 1; d < r->n; d++) { int w = (v + d) % r->n; if (r->present[w] && !r->ph[w]) { sv = w; break; } }
        int nEt = 0, keepOK = 0;
        if (pv < 0) { g->Et[nEt++] = (r->ik[v] + 1) * 1024 + v; keepOK = 1; }
        else if (!bad) {
            int ip = -1, is = -1;
            for (int k = 0; k < L; k++) if (g->al[k].k == r->tk[pv]) { ip = k; break; }
            for (int k = 0; k < L; k++) if (g->al[k].k == r->tk[sv]) { is = k; break; }
            int cnt = ((is - ip) % L + L) % L;
            if (pv == sv) cnt = L;
            for (int c = 0; c < cnt; c++) {
                int k = (ip + c) % L;
                ll key;
                if (k < L - 1) { ll a = g->al[k].k, b = g->al[k + 1].k; if (b - a < 2) { bad = 1; break; } key = a + (b - a) / 2; }
                else key = g->al[L - 1].k + KEY_GAP;
                g->Et[nEt++] = key * 1024 + v;
            }
            keepOK = pv == sv ? 1 : between(r->tk[pv], r->ik[v], r->tk[sv]);
        }
        if (bad) { respace(g, r); continue; }
        ll minik = LLONG_MAX;
        for (int u = 0; u < r->n; u++) if (r->present[u] && r->ik[u] < minik) minik = r->ik[u];
        int vbf = r->ik[v] < minik;
        Best b;
        const ll *EbUsed = g->Eb;
        if (m >= 2) {
            int nub = 0, ubs[16];
            for (int d = 1; d <= g->P->ubw; d++) for (int sg = -1; sg <= 1; sg += 2) {
                int w = v + sg * d;
                if (w < 0 || w >= r->n || g->pos[w] < 0) continue;
                int qs[2] = {g->pos[w], (g->pos[w] - 1 + m) % m};
                for (int z = 0; z < 2; z++) ubs[nub++] = qs[z];
            }
            /* slot keys increase with the slot index: sort/dedup indices, then build keys */
            for (int a = 1; a < nub; a++) { int x = ubs[a], c = a - 1; while (c >= 0 && ubs[c] > x) { ubs[c + 1] = ubs[c]; c--; } ubs[c + 1] = x; }
            int u2 = 0;
            for (int k = 0; k < nub; k++) if (u2 == 0 || ubs[k] != ubs[u2 - 1]) ubs[u2++] = ubs[k];
            nub = u2;
            for (int k = 0; k < nub; k++) {
                int q = ubs[k];
                ll key = q < m - 1 ? r->bk[g->bl[q]] + (r->bk[g->bl[q + 1]] - r->bk[g->bl[q]]) / 2 : r->bk[g->bl[m - 1]] + KEY_GAP;
                g->Eub[k] = key * 1024 + v;
            }
            Best base = solve_prep(R, r, v, NULL, 0, g->Et, nEt, keepOK, vbf, &g->st);
            ll EkI = r->ik[v] * 1024 + v;
            slot_ranks(R, ubs, nub, g->pos);
            Best ub = solve_slots(R, g->Eub, nub, nEt, vbf, EkI, base);
            R->ubInit = ub.delta == INT_MAX ? INT_MAX : ub.delta + 1;
            int ti0 = imin(imax(g->P->t0, 0), nEt - 1);
            sweep(R, r, v, m, g->pos, ti0, vbf);
            Pair *sc = xmalloc(sizeof(Pair) * m);
            for (int q = 0; q < m; q++) { sc[q].d = R->slotBest[q]; sc[q].k = q; }
            qsort(sc, m, sizeof(Pair), cmp_pair);
            int k = imin(m, g->P->topb);
            int *sel = xmalloc(sizeof(int) * (k + 1));
            for (int q = 0; q < k; q++) sel[q] = sc[q].k;
            for (int a = 1; a < k; a++) { int x = sel[a], c = a - 1; while (c >= 0 && sel[c] > x) { sel[c + 1] = sel[c]; c--; } sel[c + 1] = x; }
            for (int z = 0; z < k; z++) {
                int q = sel[z];
                ll key;
                if (q < m - 1) { ll a = r->bk[g->bl[q]], bb = r->bk[g->bl[q + 1]]; key = a + (bb - a) / 2; }
                else key = r->bk[g->bl[m - 1]] + KEY_GAP;
                g->Eb2[z] = key * 1024 + v;
            }
            R->ubInit = INT_MAX;
            int tm = R->tmode; R->tmode = 0;
            slot_ranks(R, sel, k, g->pos);
            b = solve_slots(R, g->Eb2, k, nEt, vbf, EkI, base);
            R->tmode = tm;
            EbUsed = g->Eb2;
            if (ub.delta < b.delta) { b = ub; EbUsed = g->Eub; }
            free(sc); free(sel);
        } else b = solve(R, r, v, g->Eb, nEb, g->Et, nEt, keepOK, vbf, &g->st);
        r->present[v] = 1;
        r->isA2A[v] = 0;
        if (b.keep) { r->inK[v] = 1; r->tk[v] = r->ik[v]; }
        else if (b.a2a) {
            r->inK[v] = 0; r->isA2A[v] = 1; r->tk[v] = g->Et[b.ti] >> 10;
            rep_insert_q(r, b.i, v, 2);
        } else {
            r->inK[v] = 0;
            r->bk[v] = EbUsed[b.bi] >> 10; r->tk[v] = g->Et[b.ti] >> 10;
            rep_insert_q(r, b.i, v, 0);
            rep_insert_q(r, b.j + 1, v, 1);
        }
        rep_fix_minv(r);
        return b.delta;
    }
    fprintf(stderr, "reinsert failed\n");
    exit(4);
}

static void compute_vcost(IG *g, const Rep *r) {
    for (int v = 0; v < r->n; v++) g->vcost[v] = 0;
    initial_state(r, &g->st);
    for (int k = 0; k < r->nq; k++) g->vcost[r->q[k].v] += do_step(&g->st, r, r->q[k], NULL);
}

static int pick_destroy(IG *g, const Rep *r) {
    Rng *rng = g->rng; const Params *P = g->P;
    memset(g->mark, 0, r->n);
    int nD = 0, mode = rng_below(rng, 10 + P->dbo), Q = r->nq;
    if (mode >= 10) {  /* values adjacent in BOrder */
        int m = build_bl(g, r);
        if (m > 0) {
            int L = 1 + rng_below(rng, P->seglen), st = rng_below(rng, m);
            for (int k = 0; k < imin(L, m); k++) { int v = g->bl[(st + k) % m]; if (!g->mark[v]) { g->mark[v] = 1; g->D[nD++] = v; } }
        }
    } else if (mode < 6 && Q > 0) {
        int L = 1 + rng_below(rng, P->seglen), st = rng_below(rng, Q);
        for (int k = st; k < imin(Q, st + L); k++) { int v = r->q[k].v; if (!g->mark[v]) { g->mark[v] = 1; g->D[nD++] = v; } }
    } else if (mode < 8) {
        int w = 1 + rng_below(rng, P->vwin), x = rng_below(rng, r->n);
        for (int k = 0; k < w; k++) { int v = (x + k) % r->n; if (!g->mark[v]) { g->mark[v] = 1; g->D[nD++] = v; } }
    } else {
        int c = 1 + rng_below(rng, P->rmax);
        for (int k = 0; k < c; k++) { int v = rng_below(rng, r->n); if (!g->mark[v]) { g->mark[v] = 1; g->D[nD++] = v; } }
    }
    if (nD >= r->n - 1) nD = imax(1, r->n - 2);
    return nD;
}

static IG *g_ig;
static int cmp_ord(const void *x, const void *y) {
    int a = *(const int *)x, b = *(const int *)y;
    if (g_ig->dc[a] != g_ig->dc[b]) return g_ig->dc[a] > g_ig->dc[b] ? -1 : 1;
    return (a > b) - (a < b);
}

void ig_run(const Rep *init, const Params *P, Rng *rng, Rep *best, int *bestCost, long long *itersOut) {
    int n = init->n;
    IG g; memset(&g, 0, sizeof g);
    g.P = P; g.rng = rng; g.n = n;
    g.R.tmode = 1; g.R.t0 = P->t0; g.R.topb = P->topb; g.R.ubInit = INT_MAX;
    g.R.allowA2A = P->a2a; g.R.a2aMax = P->a2amax; g.R.irot = P->irot;
    st_alloc(&g.st, n);
    g.bl = xmalloc(sizeof(int) * n); g.pos = xmalloc(sizeof(int) * n); g.vcost = xmalloc(sizeof(int) * n);
    g.D = xmalloc(sizeof(int) * n); g.dc = xmalloc(sizeof(int) * n); g.ord = xmalloc(sizeof(int) * n); g.mark = xmalloc(n);
    g.Eb = xmalloc(sizeof(ll) * (n + 2)); g.Et = xmalloc(sizeof(ll) * (2 * n + 2));
    g.Eub = xmalloc(sizeof(ll) * 16); g.Eb2 = xmalloc(sizeof(ll) * (P->topb + 2));
    g.al = xmalloc(sizeof(struct AE) * (2 * n + 2));
    Rep cur, w; memset(&cur, 0, sizeof cur); memset(&w, 0, sizeof w);
    rep_copy(&cur, init);
    int curCost = full_cost(&cur);
    rep_copy(best, &cur); *bestCost = curCost;
    int hl = imax(1, P->lahc);
    int *hist = xmalloc(sizeof(int) * hl);
    for (int k = 0; k < hl; k++) hist[k] = curCost;
    double t0 = now_sec();
    long long it;
    for (it = 0;; it++) {
        if (P->iters >= 0 ? it >= P->iters : ((it & 7) == 0 && now_sec() - t0 > P->time)) break;
        rep_copy(&w, &cur);
        compute_vcost(&g, &w);
        int nD = pick_destroy(&g, &w);
        for (int k = 0; k < nD; k++) { g.dc[k] = g.vcost[g.D[k]]; g.ord[k] = k; }
        g_ig = &g; qsort(g.ord, nD, sizeof(int), cmp_ord);
        int qd = rng_uni(rng) < P->pqd;
        for (int k = 0; k < nD; k++) {
            int v = g.D[k];
            if (qd) { rep_remove_queries(&w, v); w.inK[v] = 1; w.isA2A[v] = 0; w.tk[v] = w.ik[v]; w.ph[v] = 1; }
            else remove_value(&w, v);
        }
        int c = full_cost(&w);
        for (int z = 0; z < nD; z++) {
            int v = g.D[g.ord[z]];
            if (qd) { w.ph[v] = 0; remove_value(&w, v); c = full_cost(&w); }
            c += reinsert(&g, &w, v);
            if (P->a2a) c = full_cost(&w);  /* DP is approximate around A2A queries: use the exact cost */
        }
        int h = hist[it % hl];
        if (c <= curCost || c <= h) {
            Rep t = cur; cur = w; w = t; curCost = c;
            if (c < *bestCost) { *bestCost = c; rep_copy(best, &cur); }
        }
        hist[it % hl] = curCost;
    }
    *itersOut = it;
    rep_free(&cur); rep_free(&w); free(hist);
    st_free(&g.st);
    free(g.bl); free(g.pos); free(g.vcost); free(g.D); free(g.dc); free(g.ord); free(g.mark);
    free(g.Eb); free(g.Et); free(g.Eub); free(g.Eb2); free(g.al);
    RI *R = &g.R;
    free(R->qi); free(R->rkI); free(R->P); free(R->G); free(R->Aarr); free(R->bFbC2); free(R->bFbS2); free(R->bC1ii);
    free(R->rb); free(R->rt); free(R->tSUF); free(R->tF2); free(R->tFS); free(R->tFC); free(R->tFinA); free(R->tFinB);
    free(R->pbuf); free(R->Dg); free(R->Dr); free(R->tb); free(R->swR); free(R->swG); free(R->swPm); free(R->swPmi); free(R->slotBest); free(R->slotI); free(R->slotJ);
}
