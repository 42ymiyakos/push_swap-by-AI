/* Representation: K set, query list, AOrder (ik, tk) and BOrder (bk) keys; decoder, emitter, converter.
 * Stacks are kept cyclically sorted by key. Encoded entry = key*1024 + value. */
#include "ps.h"

void rep_alloc(Rep *r, int n) {
    r->n = n;
    r->a0 = xcalloc(n, sizeof(int));
    r->inK = xcalloc(n, 1); r->ph = xcalloc(n, 1); r->present = xcalloc(n, 1);
    r->ik = xcalloc(n, sizeof(ll)); r->tk = xcalloc(n, sizeof(ll)); r->bk = xcalloc(n, sizeof(ll));
    r->qcap = 2 * n + 8; r->q = xcalloc(r->qcap, sizeof(Query)); r->nq = 0; r->minv = 0;
}
void rep_free(Rep *r) {
    free(r->a0); free(r->inK); free(r->ph); free(r->present);
    free(r->ik); free(r->tk); free(r->bk); free(r->q);
    memset(r, 0, sizeof(*r));
}
void rep_copy(Rep *d, const Rep *s) {
    if (d->n != s->n || !d->a0) { if (d->a0) rep_free(d); rep_alloc(d, s->n); }
    int n = s->n;
    memcpy(d->a0, s->a0, sizeof(int) * n);
    memcpy(d->inK, s->inK, n); memcpy(d->ph, s->ph, n); memcpy(d->present, s->present, n);
    memcpy(d->ik, s->ik, sizeof(ll) * n); memcpy(d->tk, s->tk, sizeof(ll) * n); memcpy(d->bk, s->bk, sizeof(ll) * n);
    memcpy(d->q, s->q, sizeof(Query) * s->nq); d->nq = s->nq; d->minv = s->minv;
}
void rep_fix_minv(Rep *r) { r->minv = 0; while (r->minv < r->n && !r->present[r->minv]) r->minv++; }
void rep_insert_q(Rep *r, int pos, int v, int t) {
    memmove(r->q + pos + 1, r->q + pos, sizeof(Query) * (r->nq - pos));
    r->q[pos].v = (int16_t)v; r->q[pos].t = (uint8_t)t; r->nq++;
}
void rep_remove_queries(Rep *r, int v) {
    int k = 0;
    for (int i = 0; i < r->nq; i++) if (r->q[i].v != v) r->q[k++] = r->q[i];
    r->nq = k;
}

void st_alloc(St *s, int n) { s->A = xmalloc(sizeof(ll) * (n + 1)); s->B = xmalloc(sizeof(ll) * (n + 1)); s->na = s->nb = s->ta = s->tb = 0; }
void st_free(St *s) { free(s->A); free(s->B); }
void st_copy(St *d, const St *s) {
    memcpy(d->A, s->A, sizeof(ll) * s->na); memcpy(d->B, s->B, sizeof(ll) * s->nb);
    d->na = s->na; d->nb = s->nb; d->ta = s->ta; d->tb = s->tb;
}

static int cmp_ll(const void *a, const void *b) { ll x = *(const ll *)a, y = *(const ll *)b; return x < y ? -1 : x > y; }

void initial_state(const Rep *r, St *s) {
    s->na = 0;
    for (int i = 0; i < r->n; i++) { int x = r->a0[i]; if (r->present[x]) s->A[s->na++] = r->ik[x] * 1024 + x; }
    qsort(s->A, s->na, sizeof(ll), cmp_ll);
    int first = -1;
    for (int i = 0; i < r->n; i++) if (r->present[r->a0[i]]) { first = r->a0[i]; break; }
    s->ta = 0;
    if (first >= 0) s->ta = lb_ll(s->A, s->na, r->ik[first] * 1024 + first);
    s->nb = 0; s->tb = 0;
}

int do_step(St *s, const Rep *r, Query qu, StepInfo *info) {
    int v = qu.v;
    ll *S = qu.t == 0 ? s->A : s->B, *D = qu.t == 0 ? s->B : s->A;
    int *ns = qu.t == 0 ? &s->na : &s->nb, *nd = qu.t == 0 ? &s->nb : &s->na;
    int *ts = qu.t == 0 ? &s->ta : &s->tb, *td = qu.t == 0 ? &s->tb : &s->ta;
    ll es = (qu.t == 0 ? r->ik[v] : r->bk[v]) * 1024 + v;
    ll ed = (qu.t == 0 ? r->bk[v] : r->tk[v]) * 1024 + v;
    int ss = *ns, sd = *nd;
    int i = lb_ll(S, ss, es);
    int fs = ((i - *ts) % ss + ss) % ss;
    int j = lb_ll(D, sd, ed);
    int fd = 0;
    if (sd > 0) { int jj = j == sd ? 0 : j; fd = ((jj - *td) % sd + sd) % sd; }
    int fa = qu.t == 0 ? fs : fd, sa = qu.t == 0 ? ss : sd, fb = qu.t == 0 ? fd : fs, sb = qu.t == 0 ? sd : ss;
    int mode;
    int c = comb_cost(fa, sa, fb, sb, &mode) + 1;
    if (info) { info->fa = fa; info->sa = sa; info->fb = fb; info->sb = sb; info->mode = mode; }
    memmove(S + i, S + i + 1, sizeof(ll) * (ss - i - 1)); (*ns)--;
    *ts = *ns == 0 ? 0 : (i == *ns ? 0 : i);
    memmove(D + j + 1, D + j, sizeof(ll) * (sd - j)); D[j] = ed; (*nd)++;
    *td = j;
    return c;
}

int final_rot(const St *s, const Rep *r, int *dir) {
    int sa = s->na;
    if (sa <= 1 || r->minv >= r->n) { if (dir) *dir = 0; return 0; }
    int i = lb_ll(s->A, sa, r->tk[r->minv] * 1024 + r->minv);
    int f = ((i - s->ta) % sa + sa) % sa;
    if (dir) *dir = f <= sa - f ? 0 : 1;
    return imin(f, sa - f);
}

int full_cost(const Rep *r) {
    St s; st_alloc(&s, r->n);
    initial_state(r, &s);
    int c = 0;
    for (int k = 0; k < r->nq; k++) c += do_step(&s, r, r->q[k], NULL);
    c += final_rot(&s, r, NULL);
    st_free(&s);
    return c;
}

static void put(OpVec *ops, int o, int k) { for (int i = 0; i < k; i++) ov_push(ops, o); }

void emit(const Rep *r, OpVec *ops) {
    St s; st_alloc(&s, r->n);
    initial_state(r, &s);
    ops->n = 0;
    for (int k = 0; k < r->nq; k++) {
        StepInfo in;
        do_step(&s, r, r->q[k], &in);
        int ra = in.fa, rra = in.fa ? in.sa - in.fa : 0, rb = in.fb, rrb = in.fb ? in.sb - in.fb : 0;
        switch (in.mode) {
        case 0: put(ops, RR, imin(ra, rb)); put(ops, RA, ra - imin(ra, rb)); put(ops, RB, rb - imin(ra, rb)); break;
        case 1: put(ops, RRR, imin(rra, rrb)); put(ops, RRA, rra - imin(rra, rrb)); put(ops, RRB, rrb - imin(rra, rrb)); break;
        case 2: put(ops, RA, ra); put(ops, RRB, rrb); break;
        case 3: put(ops, RRA, rra); put(ops, RB, rb); break;
        }
        ov_push(ops, r->q[k].t == 0 ? PB : PA);
    }
    int dir, c = final_rot(&s, r, &dir);
    put(ops, dir == 0 ? RA : RRA, c);
    st_free(&s);
}

/* Ghost-list conversion of a swap-free, single-trip command sequence. */
typedef struct { int v; char live, target; } GE;
static int nextlive(const GE *L, int m, int t, int dir) {
    for (int k = 1; k <= m; k++) { int i = ((t + dir * k) % m + m) % m; if (L[i].live) return i; }
    return t;
}
static void ge_insert(GE *L, int *m, int pos, GE e) { memmove(L + pos + 1, L + pos, sizeof(GE) * (*m - pos)); L[pos] = e; (*m)++; }

int convert(const int *a0, int n, const OpVec *ops, Rep *r) {
    GE *A = xmalloc(sizeof(GE) * (2 * n + 2)), *B = xmalloc(sizeof(GE) * (n + 2));
    int mA = 0, mB = 0, ta = 0, tb = 0, liveA = n, liveB = 0, ok = 1;
    char *trip = xcalloc(n, 1);
    for (int i = 0; i < n; i++) { A[mA].v = a0[i]; A[mA].live = 1; A[mA].target = 0; mA++; }
    rep_alloc(r, n);
    memcpy(r->a0, a0, sizeof(int) * n);
    memset(r->inK, 1, n); memset(r->present, 1, n);
    for (int k = 0; k < ops->n && ok; k++) {
        int o = ops->d[k];
        if (o == SA || o == SB || o == SS) { ok = 0; break; }
        if ((o == RA || o == RR) && liveA > 1) ta = nextlive(A, mA, ta, 1);
        if ((o == RRA || o == RRR) && liveA > 1) ta = nextlive(A, mA, ta, -1);
        if ((o == RB || o == RR) && liveB > 1) tb = nextlive(B, mB, tb, 1);
        if ((o == RRB || o == RRR) && liveB > 1) tb = nextlive(B, mB, tb, -1);
        if (o == PB) {
            if (!liveA) { ok = 0; break; }
            int v = A[ta].v;
            if (A[ta].target || trip[v]) { ok = 0; break; }
            trip[v] = 1; r->inK[v] = 0; r->q[r->nq].v = (int16_t)v; r->q[r->nq].t = 0; r->nq++;
            A[ta].live = 0; liveA--;
            if (liveA) ta = nextlive(A, mA, ta, 1);
            int pos = liveB ? tb : mB;
            GE e = {v, 1, 0}; ge_insert(B, &mB, pos, e); tb = pos; liveB++;
        }
        if (o == PA) {
            if (!liveB) { ok = 0; break; }
            int v = B[tb].v; B[tb].live = 0; liveB--;
            r->q[r->nq].v = (int16_t)v; r->q[r->nq].t = 1; r->nq++;
            if (liveB) tb = nextlive(B, mB, tb, 1);
            int pos = liveA ? ta : mA;
            GE e = {v, 1, 1}; ge_insert(A, &mA, pos, e); ta = pos; liveA++;
        }
    }
    if (liveB) ok = 0;
    if (ok) {
        for (int i = 0; i < mA; i++) {
            ll key = KEY_BASE + (ll)i * KEY_GAP;
            if (A[i].target) r->tk[A[i].v] = key; else r->ik[A[i].v] = key;
        }
        for (int v = 0; v < n; v++) if (r->inK[v]) r->tk[v] = r->ik[v];
        for (int i = 0; i < mB; i++) r->bk[B[i].v] = KEY_BASE + (ll)i * KEY_GAP;
        r->minv = 0;
    }
    free(A); free(B); free(trip);
    return ok;
}
