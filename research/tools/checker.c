/* Strict push_swap checker.
 * Usage: checker <numbers...>    (ops read from stdin)
 *        checker -f input_file   (numbers from file, whitespace separated)
 * Prints "OK <count>" or "KO <reason>" / "Error". Exit 0 iff OK.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>

static int *A, *B; static int na, nb, n; /* stacks stored top at index 0 via ring */
/* simple implementation: arrays with top at end for O(1) push/pop; rotations O(n) is too slow for big
 * inputs, so use deques via ring buffers of capacity 2n */
typedef struct { int *d; int cap, head, len; } dq;
static dq a, b;
static int get(dq *q, int i) { return q->d[(q->head + i) % q->cap]; }
static void push_top(dq *q, int v) { q->head = (q->head - 1 + q->cap) % q->cap; q->d[q->head] = v; q->len++; }
static int pop_top(dq *q) { int v = q->d[q->head]; q->head = (q->head + 1) % q->cap; q->len--; return v; }
static void push_bot(dq *q, int v) { q->d[(q->head + q->len) % q->cap] = v; q->len++; }
static int pop_bot(dq *q) { q->len--; return q->d[(q->head + q->len) % q->cap]; }
static void sw(dq *q) { if (q->len < 2) return; int x = pop_top(q), y = pop_top(q); push_top(q, x); push_top(q, y); }
static void rot(dq *q) { if (q->len < 2) return; push_bot(q, pop_top(q)); }
static void rrot(dq *q) { if (q->len < 2) return; push_top(q, pop_bot(q)); }

static int parse_int(const char *s, int *out) {
    char *end; errno = 0; long v = strtol(s, &end, 10);
    if (*s == 0 || *end != 0 || errno || v < INT_MIN || v > INT_MAX) return 0;
    *out = (int)v; return 1;
}
int main(int argc, char **argv) {
    int cap = 0, *vals = NULL;
    if (argc >= 3 && !strcmp(argv[1], "-f")) {
        FILE *f = fopen(argv[2], "r"); if (!f) { puts("Error"); return 2; }
        int x; cap = 16; vals = malloc(sizeof(int) * cap);
        while (fscanf(f, "%d", &x) == 1) { if (n == cap) { cap *= 2; vals = realloc(vals, sizeof(int) * cap); } vals[n++] = x; }
        fclose(f);
    } else {
        vals = malloc(sizeof(int) * (argc + 1));
        for (int i = 1; i < argc; i++) {
            /* support single-string argument with spaces */
            char *s = strdup(argv[i]), *tok = strtok(s, " \t");
            while (tok) { int v; if (!parse_int(tok, &v)) { puts("Error"); return 2; }
                if (n >= argc + 1) { vals = realloc(vals, sizeof(int) * (n * 2 + 2)); }
                vals[n++] = v; tok = strtok(NULL, " \t"); }
            free(s);
        }
    }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (vals[i] == vals[j]) { puts("Error"); return 2; }
    int C = 2 * n + 4;
    a.d = malloc(sizeof(int) * C); b.d = malloc(sizeof(int) * C); a.cap = b.cap = C;
    for (int i = 0; i < n; i++) push_bot(&a, vals[i]);
    char line[64]; long count = 0;
    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line); if (L && line[L - 1] == '\n') line[--L] = 0;
        if (L == 0) { puts("KO empty-line"); return 1; }
        count++;
        if (!strcmp(line, "sa")) sw(&a); else if (!strcmp(line, "sb")) sw(&b);
        else if (!strcmp(line, "ss")) { sw(&a); sw(&b); }
        else if (!strcmp(line, "pa")) { if (b.len) push_top(&a, pop_top(&b)); }
        else if (!strcmp(line, "pb")) { if (a.len) push_top(&b, pop_top(&a)); }
        else if (!strcmp(line, "ra")) rot(&a); else if (!strcmp(line, "rb")) rot(&b);
        else if (!strcmp(line, "rr")) { rot(&a); rot(&b); }
        else if (!strcmp(line, "rra")) rrot(&a); else if (!strcmp(line, "rrb")) rrot(&b);
        else if (!strcmp(line, "rrr")) { rrot(&a); rrot(&b); }
        else { printf("KO bad-op:%s\n", line); return 1; }
    }
    if (b.len) { printf("KO b-not-empty %ld\n", count); return 1; }
    for (int i = 0; i + 1 < a.len; i++) if (get(&a, i) > get(&a, i + 1)) { printf("KO unsorted %ld\n", count); return 1; }
    printf("OK %ld\n", count);
    (void)A; (void)B; (void)na; (void)nb;
    return 0;
}
