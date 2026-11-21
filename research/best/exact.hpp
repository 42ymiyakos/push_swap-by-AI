// Exact optimal solver via BFS over all (A,B) states. Practical for n <= 8.
#pragma once
#include "common.hpp"
#include <unordered_map>

namespace exact {

static inline uint64_t enc(const std::vector<int> &a, const std::vector<int> &b) {
    uint64_t c = a.size();
    int sh = 4;
    for (int x : a) { c |= (uint64_t)x << sh; sh += 4; }
    for (int x : b) { c |= (uint64_t)x << sh; sh += 4; }
    return c;
}
static inline void dec(uint64_t c, int n, std::vector<int> &a, std::vector<int> &b) {
    int la = c & 15; a.clear(); b.clear(); int sh = 4;
    for (int i = 0; i < n; i++) { int x = (c >> sh) & 15; sh += 4; if (i < la) a.push_back(x); else b.push_back(x); }
}
static inline void apply(std::vector<int> &a, std::vector<int> &b, Op o) {
    auto sw = [](std::vector<int> &s) { if (s.size() >= 2) std::swap(s[0], s[1]); };
    auto rot = [](std::vector<int> &s) { if (s.size() >= 2) std::rotate(s.begin(), s.begin() + 1, s.end()); };
    auto rrot = [](std::vector<int> &s) { if (s.size() >= 2) std::rotate(s.begin(), s.end() - 1, s.end()); };
    switch (o) {
    case SA: sw(a); break; case SB: sw(b); break; case SS: sw(a); sw(b); break;
    case PA: if (!b.empty()) { a.insert(a.begin(), b[0]); b.erase(b.begin()); } break;
    case PB: if (!a.empty()) { b.insert(b.begin(), a[0]); a.erase(a.begin()); } break;
    case RA: rot(a); break; case RB: rot(b); break; case RR: rot(a); rot(b); break;
    case RRA: rrot(a); break; case RRB: rrot(b); break; case RRR: rrot(a); rrot(b); break;
    default: break;
    }
}

// ranks 0..n-1, n <= 8 (encoding uses 4 bits each, plus 4 bits length => n<=15 fits in 64 bits)
static inline std::vector<Op> solve(const std::vector<int> &init) {
    int n = init.size();
    std::vector<int> a = init, b, goal(n);
    for (int i = 0; i < n; i++) goal[i] = i;
    uint64_t s0 = enc(a, b), g = enc(goal, {});
    if (s0 == g) return {};
    std::unordered_map<uint64_t, std::pair<uint64_t, uint8_t>> par;  // state -> (parent, op)
    par.reserve(1 << 20);
    std::vector<uint64_t> cur{s0}, nxt;
    par[s0] = {s0, 255};
    std::vector<int> ta, tb;
    while (!cur.empty()) {
        nxt.clear();
        for (uint64_t c : cur) {
            for (int o = 0; o < NOPS; o++) {
                dec(c, n, ta, tb);
                apply(ta, tb, (Op)o);
                uint64_t d = enc(ta, tb);
                if (par.count(d)) continue;
                par[d] = {c, (uint8_t)o};
                if (d == g) {
                    std::vector<Op> ops;
                    for (uint64_t x = d; x != s0; x = par[x].first) ops.push_back((Op)par[x].second);
                    std::reverse(ops.begin(), ops.end());
                    return ops;
                }
                nxt.push_back(d);
            }
        }
        std::swap(cur, nxt);
    }
    return {};
}

}  // namespace exact
