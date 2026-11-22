/* Exact distance-to-sorted table for all push_swap states of size n (n <= 11), plus regret analysis.
 * State: sequence A(top..bottom) ++ B(top..bottom) and |A|; index = rank(seq) * (n+1) + |A|.
 * The op set is closed under inversion, so a BFS from the sorted state gives every state's distance.
 *   optbfs N TABLE build                      build and save the table (uint8 per state)
 *   optbfs N TABLE analyze < pairs            each input line: "v0 v1 ... v(N-1) | op op op ..."
 *     prints per case: opt, ours, and every move whose distance does not drop by 1 (wasted), with op and context. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

enum { SA, SB, SS, PA, PB, RA, RB, RR, RRA, RRB, RRR, NOPS };
static const char *NM[NOPS] = {"sa", "sb", "ss", "pa", "pb", "ra", "rb", "rr", "rra", "rrb", "rrr"};
static int N;
static uint64_t FACT[13];

static uint64_t prank(const int *s) {
    uint64_t r = 0; int used = 0;
    for (int i = 0; i < N; i++) {
        int c = 0;
        for (int j = 0; j < s[i]; j++) if (!(used >> j & 1)) c++;
        used |= 1 << s[i];
        r += (uint64_t)c * FACT[N - 1 - i];
    }
    return r;
}
static void punrank(uint64_t r, int *s) {
    int used = 0;
    for (int i = 0; i < N; i++) {
        uint64_t f = FACT[N - 1 - i]; int c = (int)(r / f); r %= f;
        for (int j = 0; j < N; j++) if (!(used >> j & 1)) { if (c-- == 0) { s[i] = j; used |= 1 << j; break; } }
    }
}
static void sw(int *s, int n) { if (n >= 2) { int t = s[0]; s[0] = s[1]; s[1] = t; } }
static void rot(int *s, int n) { if (n >= 2) { int t = s[0]; memmove(s, s + 1, sizeof(int) * (n - 1)); s[n - 1] = t; } }
static void rrot(int *s, int n) { if (n >= 2) { int t = s[n - 1]; memmove(s + 1, s, sizeof(int) * (n - 1)); s[0] = t; } }
/* a = seq[0..na-1], b = seq[na..N-1] */
static void apply(int *q, int *na, int o) {
    int a[16], b[16], la = *na, lb = N - la;
    memcpy(a, q, sizeof(int) * la); memcpy(b, q + la, sizeof(int) * lb);
    switch (o) {
    case SA: sw(a, la); break; case SB: sw(b, lb); break; case SS: sw(a, la); sw(b, lb); break;
    case PA: if (lb) { memmove(a + 1, a, sizeof(int) * la); a[0] = b[0]; memmove(b, b + 1, sizeof(int) * (lb - 1)); la++; lb--; } break;
    case PB: if (la) { memmove(b + 1, b, sizeof(int) * lb); b[0] = a[0]; memmove(a, a + 1, sizeof(int) * (la - 1)); la--; lb++; } break;
    case RA: rot(a, la); break; case RB: rot(b, lb); break; case RR: rot(a, la); rot(b, lb); break;
    case RRA: rrot(a, la); break; case RRB: rrot(b, lb); break; case RRR: rrot(a, la); rrot(b, lb); break;
    }
    memcpy(q, a, sizeof(int) * la); memcpy(q + la, b, sizeof(int) * lb); *na = la;
}
static uint64_t idx(const int *q, int na) { return prank(q) * (uint64_t)(N + 1) + (uint64_t)na; }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: optbfs N TABLE build|analyze\n"); return 1; }
    N = atoi(argv[1]);
    FACT[0] = 1; for (int i = 1; i <= 12; i++) FACT[i] = FACT[i - 1] * i;
    uint64_t S = FACT[N] * (uint64_t)(N + 1);
    uint8_t *D = malloc(S);
    if (!strcmp(argv[3], "build")) {
        memset(D, 255, S);
        uint32_t *cur = malloc(sizeof(uint32_t) * S / 4 + 16), *nxt = malloc(sizeof(uint32_t) * S / 2 + 16);
        int q[16]; for (int i = 0; i < N; i++) q[i] = i;
        uint64_t g = idx(q, N); D[g] = 0; cur[0] = (uint32_t)g; uint64_t nc = 1, tot = 1;
        for (int d = 0; nc; d++) {
            uint64_t nn = 0;
            for (uint64_t c = 0; c < nc; c++) {
                uint64_t x = cur[c]; int na = (int)(x % (N + 1)), s[16];
                punrank(x / (N + 1), s);
                for (int o = 0; o < NOPS; o++) {
                    int t[16], ta = na; memcpy(t, s, sizeof(int) * N);
                    apply(t, &ta, o);
                    uint64_t y = idx(t, ta);
                    if (D[y] == 255) { D[y] = (uint8_t)(d + 1); nxt[nn++] = (uint32_t)y; }
                }
            }
            tot += nn;
            fprintf(stderr, "depth %d: %llu states (total %llu)\n", d + 1, (unsigned long long)nn, (unsigned long long)tot);
            uint32_t *tm = cur; cur = nxt; nxt = tm; nc = nn;
            if (nc > S / 4) { fprintf(stderr, "frontier too large\n"); return 1; }
        }
        FILE *f = fopen(argv[2], "wb"); fwrite(D, 1, S, f); fclose(f);
        return 0;
    }
    FILE *f = fopen(argv[2], "rb"); if (!f || fread(D, 1, S, f) != S) { fprintf(stderr, "bad table\n"); return 1; }
    fclose(f);
    char line[1 << 16];
    long long sumOpt = 0, sumOurs = 0, cases = 0, waste[NOPS] = {0}, cnt[NOPS] = {0};
    while (fgets(line, sizeof line, stdin)) {
        char *bar = strchr(line, '|'); if (!bar) continue; *bar = 0;
        int q[16], na = N, k = 0; char *p = line, *e;
        for (; k < N; k++) { q[k] = (int)strtol(p, &e, 10); p = e; }
        int opt = D[idx(q, na)], ours = 0;
        printf("case opt %d", opt);
        char *tok = strtok(bar + 1, " \t\n");
        int prev = opt;
        char buf[1 << 14]; int bl = 0; buf[0] = 0;
        while (tok) {
            int o = -1; for (int j = 0; j < NOPS; j++) if (!strcmp(tok, NM[j])) o = j;
            if (o < 0) { tok = strtok(NULL, " \t\n"); continue; }
            apply(q, &na, o); ours++;
            int d = D[idx(q, na)];
            cnt[o]++;
            if (d != prev - 1) {
                waste[o] += d - (prev - 1);
                bl += snprintf(buf + bl, sizeof buf - bl, " %d:%s(+%d,|A|=%d)", ours, NM[o], d - (prev - 1), na);
            }
            prev = d;
            tok = strtok(NULL, " \t\n");
        }
        printf(" ours %d end %d waste%s\n", ours, prev, buf);
        sumOpt += opt; sumOurs += ours; cases++;
    }
    printf("cases %lld mean opt %.3f mean ours %.3f\n", cases, (double)sumOpt / cases, (double)sumOurs / cases);
    printf("op counts / waste:");
    for (int j = 0; j < NOPS; j++) printf(" %s %lld/%lld", NM[j], cnt[j], waste[j]);
    printf("\n");
    return 0;
}
