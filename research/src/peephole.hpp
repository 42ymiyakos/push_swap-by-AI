// Peephole optimizer: replace a window ops[i, i+K) by a shorter equivalent sequence found by
// meet-in-the-middle BFS (forward from S_i, backward from S_{i+K} with inverse ops).
// Within K ops only the K elements nearest to each end of each stack can move, so a state is
// represented by the "front" and "back" parts of each stack (the middle block is untouched).
#pragma once
#include "common.hpp"
#include <unordered_map>

namespace peep {

struct LStack {
    std::vector<int16_t> f, b;  // front (top first) and back (bottom last); mid elements implicit
    int mid = 0;                // number of untouched middle elements (0 => whole stack is in f)
};
struct LState { LStack A, B; };

static inline bool popf(LStack &s, int16_t &x) {
    if (!s.f.empty()) { x = s.f.front(); s.f.erase(s.f.begin()); return true; }
    if (s.mid == 0 && !s.b.empty()) { x = s.b.front(); s.b.erase(s.b.begin()); return true; }
    return false;  // would touch the middle block: invalid in this window
}
static inline bool popb(LStack &s, int16_t &x) {
    if (!s.b.empty()) { x = s.b.back(); s.b.pop_back(); return true; }
    if (s.mid == 0 && !s.f.empty()) { x = s.f.back(); s.f.pop_back(); return true; }
    return false;
}
static inline int size(const LStack &s) { return s.f.size() + s.b.size() + s.mid; }
static inline bool sw(LStack &s) {
    if (size(s) < 2) return true;
    int16_t x, y; if (!popf(s, x)) return false; if (!popf(s, y)) return false;
    s.f.insert(s.f.begin(), x); s.f.insert(s.f.begin(), y); return true;
}
static inline bool rot(LStack &s) {
    if (size(s) < 2) return true;
    int16_t x; if (!popf(s, x)) return false; s.b.push_back(x); return true;
}
static inline bool rrot(LStack &s) {
    if (size(s) < 2) return true;
    int16_t x; if (!popb(s, x)) return false; s.f.insert(s.f.begin(), x); return true;
}
static inline bool apply(LState &st, Op o) {
    switch (o) {
    case SA: return sw(st.A); case SB: return sw(st.B); case SS: return sw(st.A) && sw(st.B);
    case PA: { if (!size(st.B)) return true; int16_t x; if (!popf(st.B, x)) return false; st.A.f.insert(st.A.f.begin(), x); return true; }
    case PB: { if (!size(st.A)) return true; int16_t x; if (!popf(st.A, x)) return false; st.B.f.insert(st.B.f.begin(), x); return true; }
    case RA: return rot(st.A); case RB: return rot(st.B); case RR: return rot(st.A) && rot(st.B);
    case RRA: return rrot(st.A); case RRB: return rrot(st.B); case RRR: return rrot(st.A) && rrot(st.B);
    default: return false;
    }
}
static const Op INV[NOPS] = {SA, SB, SS, PB, PA, RRA, RRB, RRR, RA, RB, RR};

static inline uint64_t hs(const LState &s) {
    uint64_t h = 1469598103934665603ULL;
    auto mix = [&](uint64_t x) { h ^= x + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2); };
    for (auto *st : {&s.A, &s.B}) {
        // canonical: when mid == 0 the split between f and b is arbitrary -> hash the concatenation
        mix(0xABCD); mix(st->mid);
        if (st->mid == 0) { for (auto x : st->f) mix(x); for (auto x : st->b) mix(x); }
        else { for (auto x : st->f) mix(x); mix(0x5555); for (auto x : st->b) mix(x); }
    }
    return h;
}

static inline LStack mk(const std::deque<int> &d, int L) {
    LStack s; int n = d.size();
    if (n <= 2 * L + 2) { for (int x : d) s.f.push_back(x); return s; }
    for (int i = 0; i < L; i++) s.f.push_back(d[i]);
    for (int i = n - L; i < n; i++) s.b.push_back(d[i]);
    s.mid = n - 2 * L; return s;
}

// BFS levels: map hash -> op path (as small vector)
static inline void bfs(const LState &s0, int depth, bool inverse, std::unordered_map<uint64_t, std::vector<Op>> &seen, std::vector<uint64_t> &order) {
    std::vector<std::pair<LState, std::vector<Op>>> cur{{s0, {}}}, nxt;
    seen[hs(s0)] = {}; order.push_back(hs(s0));
    for (int d = 0; d < depth; d++) {
        nxt.clear();
        for (auto &pr : cur) for (int o = 0; o < NOPS; o++) {
            LState t = pr.first;
            if (!apply(t, (Op)o)) continue;
            uint64_t h = hs(t);
            if (seen.count(h)) continue;
            std::vector<Op> p = pr.second; p.push_back((Op)o);
            seen[h] = p; order.push_back(h);
            nxt.push_back({std::move(t), std::move(p)});
        }
        std::swap(cur, nxt);
    }
    (void)inverse;
}

// One sweep over all windows of length K. Returns number of ops saved.
static inline int optimize(const std::vector<int> &init, std::vector<Op> &ops, int K) {
    int saved = 0;
    int d1 = K / 2, d2 = (K - 1) - d1;  // total depth <= K-1
    Sim s; for (int x : init) s.a.push_back(x);
    for (size_t i = 0; i + K <= ops.size(); s.apply(ops[i]), i++) {
        LState si{mk(s.a, K + 1), mk(s.b, K + 1)};
        LState sj = si; bool ok = true;
        for (int k = 0; k < K && ok; k++) ok = apply(sj, ops[i + k]);
        if (!ok) continue;
        std::unordered_map<uint64_t, std::vector<Op>> F, B;
        std::vector<uint64_t> fo, bo;
        bfs(si, d1, false, F, fo);
        bfs(sj, d2, true, B, bo);
        int best = K; std::vector<Op> rep;
        for (uint64_t h : bo) {
            const std::vector<Op> &bp = B[h];
            auto it = F.find(h);
            if (it == F.end()) continue;
            int len = it->second.size() + bp.size();
            if (len < best) {
                best = len; rep = it->second;
                for (int k = (int)bp.size() - 1; k >= 0; k--) rep.push_back(INV[bp[k]]);
            }
        }
        if (best < K) {
            std::vector<Op> cand(ops.begin(), ops.begin() + i);
            cand.insert(cand.end(), rep.begin(), rep.end());
            cand.insert(cand.end(), ops.begin() + i + K, ops.end());
            if (verify_ops(init, cand)) { saved += (int)ops.size() - (int)cand.size(); ops.swap(cand); }
            if (i + K > ops.size()) break;
        }
    }
    return saved;
}

}  // namespace peep
