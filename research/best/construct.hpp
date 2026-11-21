// Constructive heuristics producing swap-free, single-trip command sequences.
#pragma once
#include "common.hpp"

// Longest circular increasing subsequence; returns membership mask. start chooses among optimal by rng.
static inline std::vector<char> circular_lis(const std::vector<int> &a, Rng &rng, int *len_out = nullptr) {
    int n = a.size();
    std::vector<char> best(n, 0);
    int bestlen = -1, ties = 0;
    std::vector<int> tails, tailidx, prev(n);
    for (int s = 0; s < n; s++) {
        tails.clear(); tailidx.clear();
        for (int k = 0; k < n; k++) {
            int i = (s + k) % n, x = a[i];
            int p = (int)(std::lower_bound(tails.begin(), tails.end(), x) - tails.begin());
            prev[i] = p ? tailidx[p - 1] : -1;
            if (p == (int)tails.size()) { tails.push_back(x); tailidx.push_back(i); }
            else { tails[p] = x; tailidx[p] = i; }
        }
        int L = tails.size();
        bool take = false;
        if (L > bestlen) { bestlen = L; ties = 1; take = true; }
        else if (L == bestlen) { ties++; take = rng.below(ties) == 0; }
        if (take) {
            std::fill(best.begin(), best.end(), 0);
            for (int i = tailidx.back(); i >= 0; i = prev[i]) best[a[i]] = 1;
        }
    }
    if (len_out) *len_out = bestlen;
    return best;
}

// Longest increasing subsequence reading A cyclically from position s (deterministic).
static inline std::vector<char> lis_from(const std::vector<int> &a, int s) {
    int n = a.size();
    std::vector<int> tails, tailidx, prev(n);
    for (int k = 0; k < n; k++) {
        int i = (s + k) % n, x = a[i];
        int p = (int)(std::lower_bound(tails.begin(), tails.end(), x) - tails.begin());
        prev[i] = p ? tailidx[p - 1] : -1;
        if (p == (int)tails.size()) { tails.push_back(x); tailidx.push_back(i); }
        else { tails[p] = x; tailidx[p] = i; }
    }
    std::vector<char> m(n, 0);
    for (int i = tailidx.back(); i >= 0; i = prev[i]) m[a[i]] = 1;
    return m;
}

struct Builder {
    Sim s; std::vector<Op> ops;
    void op(Op o, int k = 1) { for (int i = 0; i < k; i++) { s.apply(o); ops.push_back(o); } }
};

// Chunk construction: non-K values are pushed to B in value-range chunks built middle-out
// (top chunk / bottom chunk pairs), then greedily inserted back with minimal combined rotation.
#include "beam.hpp"
static inline std::vector<Op> chunk_construct(const std::vector<int> &a0, const std::vector<char> &inK, int m, int bw = 0, int bc = 2) {
    int n = a0.size();
    Builder b; for (int x : a0) b.s.a.push_back(x);
    auto chunk = [&](int v) { return m - 1 - (int)((long long)v * m / n); };  // 0 = largest values
    // rounds: list of (top chunk, bottom chunk); -1 if none
    std::vector<std::pair<int, int>> rounds;
    if (m % 2 == 1) { rounds.push_back({m / 2, -1}); for (int r = 1; r <= m / 2; r++) rounds.push_back({m / 2 - r, m / 2 + r}); }
    else for (int r = 0; r < m / 2; r++) rounds.push_back({m / 2 - 1 - r, m / 2 + r});
    int pending = 0;
    for (auto [tc, bc] : rounds) {
        while (true) {
            int sa = b.s.a.size(), bi = -1, bd = INT_MAX, bdir = 0;
            for (int i = 0; i < sa; i++) {
                int v = b.s.a[i];
                if (inK[v]) continue;
                int c = chunk(v);
                if (c != tc && c != bc) continue;
                int top = c == tc;
                // forward cost (merge with pending rb)
                int f = i, cf = top ? std::max(f, pending) : f;  // rr merges pending rb for free
                int rv = (sa - i) % sa, cr = rv + (top ? pending : 0);
                int d = std::min(cf, cr);
                if (d < bd) { bd = d; bi = i; bdir = cf <= cr ? 0 : 1; }
            }
            if (bi < 0) break;
            int v = b.s.a[bi];
            bool top = chunk(v) == tc;
            if (bdir == 0) {
                int f = bi, k = std::min(f, pending);
                b.op(RR, k); pending -= k; b.op(RA, f - k);
                if (top) { b.op(RB, pending); pending = 0; }
            } else {
                b.op(RRA, (sa - bi) % sa);
                if (top) { b.op(RB, pending); pending = 0; }
            }
            b.op(PB);
            if (!top) pending++;
        }
    }
    if (bw > 0 && n <= 512) {
        std::vector<int> Av(b.s.a.begin(), b.s.a.end()), Bv(b.s.b.begin(), b.s.b.end());
        std::vector<Op> rest = beam::solve(n, Av, b.s.a.empty() ? 0 : b.s.a[0], Bv, bw, bc);
        for (Op o : rest) b.op(o);
        return b.ops;
    }
    // B -> A greedy insertion
    std::vector<int> idxA(n), succ(n);
    while (!b.s.b.empty()) {
        int sa = b.s.a.size(), sb = b.s.b.size();
        std::vector<char> inA(n, 0);
        for (int i = 0; i < sa; i++) { idxA[b.s.a[i]] = i; inA[b.s.a[i]] = 1; }
        int nxt = -1;
        for (int v = n - 1; v >= 0; v--) { succ[v] = nxt; if (inA[v]) nxt = v; }
        int minA = nxt;  // smallest value in A
        int bestc = INT_MAX, bj = -1, bfa = 0;
        for (int j = 0; j < sb; j++) {
            int v = b.s.b[j];
            int s2 = succ[v] >= 0 ? succ[v] : minA;
            int fa = sa ? idxA[s2] : 0;
            int c = comb_cost_simple(fa, sa, j, sb);
            if (c < bestc) { bestc = c; bj = j; bfa = fa; }
        }
        int fa = bfa, fb = bj;
        int ra = fa, rra = fa ? sa - fa : 0, rb = fb, rrb = fb ? sb - fb : 0;
        int c0 = std::max(ra, rb), c1 = std::max(rra, rrb), c2 = ra + rrb, c3 = rra + rb;
        int mc = std::min(std::min(c0, c1), std::min(c2, c3));
        if (mc == c0) { int k = std::min(ra, rb); b.op(RR, k); b.op(RA, ra - k); b.op(RB, rb - k); }
        else if (mc == c1) { int k = std::min(rra, rrb); b.op(RRR, k); b.op(RRA, rra - k); b.op(RRB, rrb - k); }
        else if (mc == c2) { b.op(RA, ra); b.op(RRB, rrb); }
        else { b.op(RRA, rra); b.op(RB, rb); }
        b.op(PA);
    }
    int sa = b.s.a.size(), i0 = 0;
    for (int i = 0; i < sa; i++) if (b.s.a[i] == 0) i0 = i;
    if (i0 <= sa - i0) b.op(RA, i0); else b.op(RRA, sa - i0);
    return b.ops;
}
