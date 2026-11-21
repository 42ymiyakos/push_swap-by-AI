/* push_swap solver (C port of PSIG: query representation + Iterated Greedy with exact DP reinsertion).
 * Not norminette-compliant by design. Build: make -C research/c
 */
#ifndef PS_H
#define PS_H

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;

enum { SA, SB, SS, PA, PB, RA, RB, RR, RRA, RRB, RRR, NOPS };
extern const char *OPNAME[NOPS];

/* ---------- small utilities ---------- */
void *xmalloc(size_t n);
void *xcalloc(size_t n, size_t sz);
void *xrealloc(void *p, size_t n);
static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

/* growable op list */
typedef struct { uint8_t *d; int n, cap; } OpVec;
void ov_init(OpVec *v);
void ov_free(OpVec *v);
void ov_push(OpVec *v, int o);
void ov_copy(OpVec *dst, const OpVec *src);

/* rng (xorshift, identical to the C++ version) */
typedef struct { uint64_t s; } Rng;
void rng_seed(Rng *r, uint64_t seed);
uint64_t rng_next(Rng *r);
int rng_below(Rng *r, int n);
double rng_uni(Rng *r);

/* exact simulator (ring buffers) used for verification and constructors */
typedef struct { int *d; int cap, head, len; } Dq;
typedef struct { Dq a, b; } Sim;
void sim_init(Sim *s, const int *init, int n);
void sim_free(Sim *s);
void sim_copy(Sim *dst, const Sim *src);
void sim_apply(Sim *s, int op);
static inline int dq_get(const Dq *q, int i) { return q->d[(q->head + i) % q->cap]; }
int sim_sorted(const Sim *s);
int verify_ops(const int *init, int n, const OpVec *ops);

int lb_ll(const ll *a, int n, ll key);  /* lower_bound */
int lb_int(const int *a, int n, int key);
double now_sec(void);
void print_ops(const OpVec *ops);

static inline int comb_cost(int fa, int sa, int fb, int sb, int *mode) {
    int ra = fa, rra = fa ? sa - fa : 0, rb = fb, rrb = fb ? sb - fb : 0;
    int c0 = imax(ra, rb), c1 = imax(rra, rrb), c2 = ra + rrb, c3 = rra + rb;
    int best = c0, m = 0;
    if (c1 < best) { best = c1; m = 1; }
    if (c2 < best) { best = c2; m = 2; }
    if (c3 < best) { best = c3; m = 3; }
    if (mode) *mode = m;
    return best;
}

/* push cost with pending A rotations (pr = ra, prr = rra) of a preceding A2A run that can be paired
 * with this push's B rotations (rr / rrr). */
static inline int comb_sync(int fa, int sa, int fb, int sb, int pr, int prr, int *mode, int *kk) {
    int ra = fa, rra = fa ? sa - fa : 0, rb = fb, rrb = fb ? sb - fb : 0;
    int k0 = imin(pr, imax(0, rb - ra)), k1 = imin(prr, imax(0, rrb - rra));
    int k2 = imin(prr, rrb), k3 = imin(pr, rb);
    int c0 = imax(ra, rb) - k0, c1 = imax(rra, rrb) - k1, c2 = ra + rrb - k2, c3 = rra + rb - k3;
    int best = c0, m = 0, k = k0;
    if (c1 < best) { best = c1; m = 1; k = k1; }
    if (c2 < best) { best = c2; m = 2; k = k2; }
    if (c3 < best) { best = c3; m = 3; k = k3; }
    if (mode) *mode = m;
    if (kk) *kk = k;
    return best;
}

/* ---------- representation ---------- */
#define KEY_BASE (1LL << 40)
#define KEY_GAP (1LL << 22)

typedef struct { int16_t v; uint8_t t; } Query; /* t: 0 = A2B, 1 = B2A, 2 = A2A (move inside A with sa) */

typedef struct {
    int n;
    int *a0;
    char *inK, *ph, *present, *isA2A;
    ll *ik, *tk, *bk;
    Query *q;
    int nq, qcap;
    int minv;
} Rep;

void rep_alloc(Rep *r, int n);
void rep_free(Rep *r);
void rep_copy(Rep *dst, const Rep *src);
void rep_fix_minv(Rep *r);
void rep_insert_q(Rep *r, int pos, int v, int t);
void rep_remove_queries(Rep *r, int v);

typedef struct { ll *A, *B; int na, nb, ta, tb, pr, prr; } St;
void st_alloc(St *s, int n);
void st_free(St *s);
void st_copy(St *dst, const St *src);
void initial_state(const Rep *r, St *s);
typedef struct { int fa, sa, fb, sb, mode, k; } StepInfo;
typedef struct { int cost, fwd, steps, rot, rotfwd, topWo; } A2AInfo;
A2AInfo a2a_cost(int sa, int ta, int src, int dstRaw);
int do_step_a2a(St *s, const Rep *r, Query qu, A2AInfo *info);
int do_step(St *s, const Rep *r, Query qu, StepInfo *info);
int final_rot(const St *s, const Rep *r, int *dir);
int full_cost(const Rep *r);
void emit(const Rep *r, OpVec *ops);
int convert(const int *a0, int n, const OpVec *ops, Rep *r);
static inline int topv(const ll *s, int n, int t) { return n == 0 ? -1 : (int)(s[t] & 1023); }

/* ---------- construction ---------- */
void lis_from(const int *a, int n, int s, char *mask);
void circular_lis(const int *a, int n, Rng *rng, char *mask);
extern int g_cmode;
void chunk_construct(const int *a0, int n, const char *inK, int m, int bw, int bc, OpVec *out);
void beam_solve(int n, const int *Avals, int na, int topA, const int *Bvals, int nbv, int width, int cand, OpVec *out);

/* ---------- exact / peephole ---------- */
void exact_solve(const int *init, int n, OpVec *out);
int peep_optimize(const int *init, int n, OpVec *ops, int K);

/* ---------- iterated greedy ---------- */
typedef struct {
    double time;
    long long iters;
    uint64_t seed;
    int t0, topb, seglen, rmax, vwin, lahc, verbose;
    double pqd;
    int a2a, a2amax, ubw, irot, dbo, dbias, kick, kickmul, kickbest;
} Params;

typedef struct { int delta, i, j, bi, ti, keep, a2a; } Best;

void ig_run(const Rep *init, const Params *P, Rng *rng, Rep *best, int *bestCost, long long *itersOut);

#endif
