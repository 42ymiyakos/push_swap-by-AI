/* Exact optimum by BFS over all (A,B) states, n <= 8 (4-bit encoding). */
#include "ps.h"

static uint64_t enc(const int *a, int na, const int *b, int nb) {
    uint64_t c = (uint64_t)na; int sh = 4;
    for (int i = 0; i < na; i++) { c |= (uint64_t)a[i] << sh; sh += 4; }
    for (int i = 0; i < nb; i++) { c |= (uint64_t)b[i] << sh; sh += 4; }
    return c;
}
static void dec(uint64_t c, int n, int *a, int *na, int *b, int *nb) {
    int la = (int)(c & 15), sh = 4; *na = 0; *nb = 0;
    for (int i = 0; i < n; i++) { int x = (int)((c >> sh) & 15); sh += 4; if (i < la) a[(*na)++] = x; else b[(*nb)++] = x; }
}
static void sw(int *s, int n) { if (n >= 2) { int t = s[0]; s[0] = s[1]; s[1] = t; } }
static void rot(int *s, int n) { if (n >= 2) { int t = s[0]; memmove(s, s + 1, sizeof(int) * (n - 1)); s[n - 1] = t; } }
static void rrot(int *s, int n) { if (n >= 2) { int t = s[n - 1]; memmove(s + 1, s, sizeof(int) * (n - 1)); s[0] = t; } }
static void push(int *from, int *nf, int *to, int *nt) {
    if (!*nf) return;
    int x = from[0]; memmove(from, from + 1, sizeof(int) * (*nf - 1)); (*nf)--;
    memmove(to + 1, to, sizeof(int) * *nt); to[0] = x; (*nt)++;
}
static void apply(int *a, int *na, int *b, int *nb, int o) {
    switch (o) {
    case SA: sw(a, *na); break; case SB: sw(b, *nb); break; case SS: sw(a, *na); sw(b, *nb); break;
    case PA: push(b, nb, a, na); break; case PB: push(a, na, b, nb); break;
    case RA: rot(a, *na); break; case RB: rot(b, *nb); break; case RR: rot(a, *na); rot(b, *nb); break;
    case RRA: rrot(a, *na); break; case RRB: rrot(b, *nb); break; case RRR: rrot(a, *na); rrot(b, *nb); break;
    }
}

typedef struct { uint64_t key, par; uint8_t op, used; } Slot;

void exact_solve(const int *init, int n, OpVec *out) {
    out->n = 0;
    int a[16], b[16], na, nb, goal[16];
    for (int i = 0; i < n; i++) goal[i] = i;
    uint64_t s0 = enc(init, n, NULL, 0), g = enc(goal, n, NULL, 0);
    if (s0 == g) return;
    size_t cap = 1 << 20, mask = cap - 1;
    Slot *T = xcalloc(cap, sizeof(Slot));
    uint64_t *cur = xmalloc(sizeof(uint64_t) * 400000), *nxt = xmalloc(sizeof(uint64_t) * 400000);
    size_t ncur = 1, nnxt;
    cur[0] = s0;
#define FIND(k) ({ size_t h_ = (size_t)(((k) * 0x9E3779B97F4A7C15ULL) >> 40) & mask; while (T[h_].used && T[h_].key != (k)) h_ = (h_ + 1) & mask; h_; })
    { size_t h = FIND(s0); T[h].used = 1; T[h].key = s0; T[h].par = s0; T[h].op = 255; }
    while (ncur) {
        nnxt = 0;
        for (size_t ci = 0; ci < ncur; ci++) {
            uint64_t c = cur[ci];
            for (int o = 0; o < NOPS; o++) {
                dec(c, n, a, &na, b, &nb);
                apply(a, &na, b, &nb, o);
                uint64_t d = enc(a, na, b, nb);
                size_t h = FIND(d);
                if (T[h].used) continue;
                T[h].used = 1; T[h].key = d; T[h].par = c; T[h].op = (uint8_t)o;
                if (d == g) {
                    for (uint64_t x = d; x != s0;) { size_t hx = FIND(x); ov_push(out, T[hx].op); x = T[hx].par; }
                    for (int i = 0, j = out->n - 1; i < j; i++, j--) { uint8_t t = out->d[i]; out->d[i] = out->d[j]; out->d[j] = t; }
                    free(T); free(cur); free(nxt);
                    return;
                }
                nxt[nnxt++] = d;
            }
        }
        uint64_t *t = cur; cur = nxt; nxt = t; ncur = nnxt;
    }
#undef FIND
    free(T); free(cur); free(nxt);
}
