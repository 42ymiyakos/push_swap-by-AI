// Beam search for the B->A phase when A holds only settled (cyclically sorted) values.
// State: bitset of values in A, bitset of B positions still present, top value of A, top position in B.
#pragma once
#include "common.hpp"
#include <array>

namespace beam {

static const int W = 8;  // 64-bit words, n <= 512
typedef std::array<uint64_t, W> BS;

static inline void setb(BS &b, int i) { b[i >> 6] |= 1ULL << (i & 63); }
static inline void clrb(BS &b, int i) { b[i >> 6] &= ~(1ULL << (i & 63)); }
static inline bool getb(const BS &b, int i) { return b[i >> 6] >> (i & 63) & 1; }
// count of set bits in [l, r)
static inline int cnt(const BS &b, int l, int r) {
    if (l >= r) return 0;
    int c = 0, wl = l >> 6, wr = (r - 1) >> 6;
    if (wl == wr) return __builtin_popcountll(b[wl] >> (l & 63) & (~0ULL >> (63 - ((r - 1) & 63) + (l & 63))));
    c += __builtin_popcountll(b[wl] >> (l & 63));
    for (int w = wl + 1; w < wr; w++) c += __builtin_popcountll(b[w]);
    c += __builtin_popcountll(b[wr] & (~0ULL >> (63 - ((r - 1) & 63))));
    return c;
}
// cyclic count in [l, r) over [0, n)
static inline int ccnt(const BS &b, int l, int r, int n) { return l <= r ? cnt(b, l, r) : cnt(b, l, n) + cnt(b, 0, r); }
static inline int next_set(const BS &b, int i, int n) {  // first set index >= i cyclic, -1 if none
    for (int k = 0; k < 2; k++) {
        int lo = k == 0 ? i : 0, hi = k == 0 ? n : i;
        for (int x = lo; x < hi;) {
            int w = x >> 6; uint64_t m = b[w] >> (x & 63);
            if (m) { int y = x + __builtin_ctzll(m); return y < hi ? y : -1; }
            x = (w + 1) << 6;
        }
    }
    return -1;
}
static inline int prev_set(const BS &b, int i, int n) {  // last set index <= i cyclic
    for (int k = 0; k < n; k++) { int x = ((i - k) % n + n) % n; if (getb(b, x)) return x; }
    return -1;
}

struct Node { BS a, b; int ta, tb, cost, nA, nB; int parent; int16_t v; uint8_t mode; int fa, fb; float ev; };

// Bvals: B contents top-first. A values: settled (sorted cyclically) with given top value.
// Returns ops for the B->A phase including final rotation.
static inline std::vector<Op> solve(int n, const std::vector<int> &Avals, int topA, const std::vector<int> &Bvals, int width,
                                    int cand) {
    int m = Bvals.size();
    Node root{};
    for (int x : Avals) setb(root.a, x);
    for (int i = 0; i < m; i++) setb(root.b, i);
    root.ta = topA; root.tb = 0; root.cost = 0; root.nA = Avals.size(); root.nB = m; root.parent = -1;
    std::vector<Node> all; all.reserve((size_t)width * cand * (m + 1) + 16);
    all.push_back(root);
    std::vector<int> layer{0}, nxt;
    struct C { float ev; int node; };
    // bitset of B values still in B is derived per state from positions
    auto heur = [&](const Node &s) {
        int L = s.nA; if (L <= 1) return 0.0f;
        BS bv{};
        for (int i = 0; i < m; i++) if (getb(s.b, i)) setb(bv, Bvals[i]);
        // marked ranks: A elements whose cyclic predecessor gap contains a B value
        std::vector<int> pts; pts.reserve(64);
        int first = next_set(s.a, 0, n), prevv = -1, rk = 0, st = 0, e = 0;
        int lastv = prev_set(s.a, n - 1, n);
        for (int x = first; x >= 0;) {
            if (x == s.ta) st = rk;
            int lo = prevv < 0 ? lastv : prevv;  // gap (lo, x) cyclic
            bool mk = lo < x ? cnt(bv, lo + 1, x) > 0 : (cnt(bv, lo + 1, n) + cnt(bv, 0, x)) > 0;
            if (mk) pts.push_back(rk);
            prevv = x; rk++;
            int nx = x + 1 < n ? next_set(s.a, x + 1, n) : -1;
            if (nx <= x) break;
            x = nx;
        }
        if (pts.empty()) { int d = st; return (float)std::min(d, L - d); }
        // final target: rank 0 (minimum), unless min value in B -> its slot = rank 0 as well
        pts.push_back(st); pts.push_back(e);
        std::sort(pts.begin(), pts.end());
        pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
        int k = pts.size();
        int gmax = -1, gi = 0;
        for (int i = 0; i < k; i++) {
            int g = (i + 1 < k ? pts[i + 1] : pts[0] + L) - pts[i];
            if (g > gmax) { gmax = g; gi = i; }
        }
        int l = pts[(gi + 1) % k], arc = L - gmax;
        auto off = [&](int x) { return ((x - l) % L + L) % L; };
        int os = off(st), oe = off(e);
        int c1 = os + arc + (arc - oe), c2 = (arc - os) + arc + oe;
        return (float)std::min(c1, c2);
    };
    for (int depth = 0; depth < m; depth++) {
        std::vector<C> cands;
        for (int id : layer) {
            const Node s = all[id];
            int sa = s.nA, sb = s.nB;
            // enumerate present B elements by distance from top
            struct E { int c, pos, fa, fb, mode; };
            std::vector<E> es;
            int pf = s.tb, pb = s.tb;
            int lim = std::min(sb, 2 * cand + 40);
            std::vector<int> seen;
            // forward walk
            int dpos = 0;
            for (int x = pf, k = 0; k < sb && k <= lim; k++) {
                x = next_set(s.b, x, m); if (x < 0) break;
                seen.push_back(x); x = (x + 1) % m; dpos++;
            }
            for (int x = pb, k = 0; k < sb && k <= lim; k++) {
                x = prev_set(s.b, (x - 1 + m) % m, m); if (x < 0) break;
                seen.push_back(x);
            }
            std::sort(seen.begin(), seen.end()); seen.erase(std::unique(seen.begin(), seen.end()), seen.end());
            for (int pos : seen) {
                int v = Bvals[pos];
                int fb = ccnt(s.b, s.tb, pos, m);
                int sv = next_set(s.a, v, n);
                int fa = ccnt(s.a, s.ta, sv, n);
                int ra = fa, rra = fa ? sa - fa : 0, rb = fb, rrb = fb ? sb - fb : 0;
                int c0 = std::max(ra, rb), c1 = std::max(rra, rrb), c2 = ra + rrb, c3 = rra + rb;
                int c = c0, md = 0;
                if (c1 < c) c = c1, md = 1; if (c2 < c) c = c2, md = 2; if (c3 < c) c = c3, md = 3;
                es.push_back({c, pos, fa, fb, md});
            }
            std::sort(es.begin(), es.end(), [](const E &x, const E &y) { return x.c < y.c; });
            int take = 0;
            for (size_t k = 0; k < es.size(); k++) {
                if (take >= cand && es[k].c > es[take - 1].c) break;
                take++;
                const E &e = es[k];
                Node ch = s;
                int v = Bvals[e.pos];
                clrb(ch.b, e.pos); ch.nB--;
                setb(ch.a, v); ch.nA++; ch.ta = v;
                ch.tb = ch.nB ? next_set(ch.b, (e.pos + 1) % m, m) : 0;
                ch.cost = s.cost + e.c + 1; ch.parent = id; ch.v = v; ch.mode = e.mode; ch.fa = e.fa; ch.fb = e.fb;
                ch.ev = ch.cost + heur(ch);
                all.push_back(ch);
                cands.push_back({ch.ev, (int)all.size() - 1});
            }
        }
        // dedup by (a,b,ta,tb) hash
        std::sort(cands.begin(), cands.end(), [](const C &x, const C &y) { return x.ev < y.ev; });
        nxt.clear();
        std::vector<uint64_t> hs;
        for (auto &c : cands) {
            const Node &s = all[c.node];
            uint64_t h = 1469598103934665603ULL;
            for (int w = 0; w < W; w++) { h ^= s.b[w]; h *= 1099511628211ULL; }
            h ^= (uint64_t)s.ta * 1315423911u + s.tb; h *= 1099511628211ULL;
            if (std::find(hs.begin(), hs.end(), h) != hs.end()) continue;
            hs.push_back(h); nxt.push_back(c.node);
            if ((int)nxt.size() >= width) break;
        }
        layer.swap(nxt);
    }
    // pick best final (cost + final rotation)
    int bestId = -1, bestC = INT_MAX, bestF = 0;
    for (int id : layer) {
        const Node &s = all[id];
        int f = cnt(s.a, 0, s.ta) /* rank of top */;
        int L = s.nA, fr = std::min(f, (L - f) % L);
        if (s.cost + fr < bestC) { bestC = s.cost + fr; bestId = id; bestF = f; }
    }
    std::vector<int> path;
    for (int id = bestId; id > 0; id = all[id].parent) path.push_back(id);
    std::reverse(path.begin(), path.end());
    std::vector<Op> ops;
    auto put = [&](Op o, int k) { for (int i = 0; i < k; i++) ops.push_back(o); };
    for (int id : path) {
        const Node &c = all[id], &p = all[c.parent];
        int sa = p.nA, sb = p.nB;
        int ra = c.fa, rra = c.fa ? sa - c.fa : 0, rb = c.fb, rrb = c.fb ? sb - c.fb : 0;
        switch (c.mode) {
        case 0: put(RR, std::min(ra, rb)); put(RA, ra - std::min(ra, rb)); put(RB, rb - std::min(ra, rb)); break;
        case 1: put(RRR, std::min(rra, rrb)); put(RRA, rra - std::min(rra, rrb)); put(RRB, rrb - std::min(rra, rrb)); break;
        case 2: put(RA, ra); put(RRB, rrb); break;
        case 3: put(RRA, rra); put(RB, rb); break;
        }
        ops.push_back(PA);
    }
    const Node &s = all[bestId];
    int L = s.nA;
    if (bestF <= L - bestF) put(RRA, bestF); else put(RA, L - bestF);
    return ops;
}

}  // namespace beam
