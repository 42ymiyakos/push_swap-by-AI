#include "ps.h"
#include <time.h>

const char *OPNAME[NOPS] = {"sa", "sb", "ss", "pa", "pb", "ra", "rb", "rr", "rra", "rrb", "rrr"};

void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) { fprintf(stderr, "Error\n"); exit(1); }
    return p;
}
void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n ? n : 1, sz ? sz : 1);
    if (!p) { fprintf(stderr, "Error\n"); exit(1); }
    return p;
}
void *xrealloc(void *p, size_t n) {
    p = realloc(p, n ? n : 1);
    if (!p) { fprintf(stderr, "Error\n"); exit(1); }
    return p;
}

void ov_init(OpVec *v) { v->d = NULL; v->n = v->cap = 0; }
void ov_free(OpVec *v) { free(v->d); ov_init(v); }
void ov_push(OpVec *v, int o) {
    if (v->n == v->cap) { v->cap = v->cap ? v->cap * 2 : 64; v->d = xrealloc(v->d, v->cap); }
    v->d[v->n++] = (uint8_t)o;
}
void ov_copy(OpVec *dst, const OpVec *src) {
    if (dst->cap < src->n) { dst->cap = src->n; dst->d = xrealloc(dst->d, dst->cap); }
    if (src->n) memcpy(dst->d, src->d, src->n);
    dst->n = src->n;
}

void rng_seed(Rng *r, uint64_t seed) { r->s = seed ? seed : 88172645463325252ULL; }
uint64_t rng_next(Rng *r) { r->s ^= r->s << 7; r->s ^= r->s >> 9; return r->s; }
int rng_below(Rng *r, int n) { return (int)((rng_next(r) >> 11) % (uint64_t)n); }
double rng_uni(Rng *r) { return (rng_next(r) >> 11) * (1.0 / 9007199254740992.0); }

/* ring-buffer deque */
static void dq_init(Dq *q, int cap) { q->cap = cap + 1; q->d = xmalloc(sizeof(int) * q->cap); q->head = 0; q->len = 0; }
static void dq_push_top(Dq *q, int v) { q->head = (q->head - 1 + q->cap) % q->cap; q->d[q->head] = v; q->len++; }
static int dq_pop_top(Dq *q) { int v = q->d[q->head]; q->head = (q->head + 1) % q->cap; q->len--; return v; }
static void dq_push_bot(Dq *q, int v) { q->d[(q->head + q->len) % q->cap] = v; q->len++; }
static int dq_pop_bot(Dq *q) { q->len--; return q->d[(q->head + q->len) % q->cap]; }
static void dq_sw(Dq *q) {
    if (q->len < 2) return;
    int i0 = q->head, i1 = (q->head + 1) % q->cap, t = q->d[i0];
    q->d[i0] = q->d[i1]; q->d[i1] = t;
}
static void dq_rot(Dq *q) { if (q->len >= 2) dq_push_bot(q, dq_pop_top(q)); }
static void dq_rrot(Dq *q) { if (q->len >= 2) dq_push_top(q, dq_pop_bot(q)); }

void sim_init(Sim *s, const int *init, int n) {
    dq_init(&s->a, n); dq_init(&s->b, n);
    for (int i = 0; i < n; i++) dq_push_bot(&s->a, init[i]);
}
void sim_free(Sim *s) { free(s->a.d); free(s->b.d); }
void sim_copy(Sim *dst, const Sim *src) {
    dst->a = src->a; dst->b = src->b;
    dst->a.d = xmalloc(sizeof(int) * src->a.cap); memcpy(dst->a.d, src->a.d, sizeof(int) * src->a.cap);
    dst->b.d = xmalloc(sizeof(int) * src->b.cap); memcpy(dst->b.d, src->b.d, sizeof(int) * src->b.cap);
}
void sim_apply(Sim *s, int o) {
    switch (o) {
    case SA: dq_sw(&s->a); break;
    case SB: dq_sw(&s->b); break;
    case SS: dq_sw(&s->a); dq_sw(&s->b); break;
    case PA: if (s->b.len) dq_push_top(&s->a, dq_pop_top(&s->b)); break;
    case PB: if (s->a.len) dq_push_top(&s->b, dq_pop_top(&s->a)); break;
    case RA: dq_rot(&s->a); break;
    case RB: dq_rot(&s->b); break;
    case RR: dq_rot(&s->a); dq_rot(&s->b); break;
    case RRA: dq_rrot(&s->a); break;
    case RRB: dq_rrot(&s->b); break;
    case RRR: dq_rrot(&s->a); dq_rrot(&s->b); break;
    default: break;
    }
}
int sim_sorted(const Sim *s) {
    if (s->b.len) return 0;
    for (int i = 1; i < s->a.len; i++) if (dq_get(&s->a, i - 1) > dq_get(&s->a, i)) return 0;
    return 1;
}
int verify_ops(const int *init, int n, const OpVec *ops) {
    Sim s; sim_init(&s, init, n);
    for (int i = 0; i < ops->n; i++) sim_apply(&s, ops->d[i]);
    int ok = sim_sorted(&s);
    sim_free(&s);
    return ok;
}

int lb_ll(const ll *a, int n, ll key) {
    int lo = 0, hi = n;
    while (lo < hi) { int mid = (lo + hi) >> 1; if (a[mid] < key) lo = mid + 1; else hi = mid; }
    return lo;
}
int lb_int(const int *a, int n, int key) {
    int lo = 0, hi = n;
    while (lo < hi) { int mid = (lo + hi) >> 1; if (a[mid] < key) lo = mid + 1; else hi = mid; }
    return lo;
}
double now_sec(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}
void print_ops(const OpVec *ops) {
    size_t cap = (size_t)ops->n * 4 + 1, len = 0;
    char *buf = xmalloc(cap);
    for (int i = 0; i < ops->n; i++) {
        const char *s = OPNAME[ops->d[i]];
        while (*s) buf[len++] = *s++;
        buf[len++] = '\n';
    }
    fwrite(buf, 1, len, stdout);
    free(buf);
}
