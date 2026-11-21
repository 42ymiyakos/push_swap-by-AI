// Common utilities: operations, exact simulator, parsing.
#pragma once
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <string>
#include <vector>

enum Op : uint8_t { SA, SB, SS, PA, PB, RA, RB, RR, RRA, RRB, RRR, NOPS };
static const char *OPNAME[NOPS] = {"sa", "sb", "ss", "pa", "pb", "ra", "rb", "rr", "rra", "rrb", "rrr"};

// Exact reference simulator (deque based). top = front.
struct Sim {
    std::deque<int> a, b;
    static void sw(std::deque<int> &s) { if (s.size() >= 2) std::swap(s[0], s[1]); }
    static void rot(std::deque<int> &s) { if (s.size() >= 2) { s.push_back(s.front()); s.pop_front(); } }
    static void rrot(std::deque<int> &s) { if (s.size() >= 2) { s.push_front(s.back()); s.pop_back(); } }
    void apply(Op o) {
        switch (o) {
        case SA: sw(a); break; case SB: sw(b); break; case SS: sw(a); sw(b); break;
        case PA: if (!b.empty()) { a.push_front(b.front()); b.pop_front(); } break;
        case PB: if (!a.empty()) { b.push_front(a.front()); a.pop_front(); } break;
        case RA: rot(a); break; case RB: rot(b); break; case RR: rot(a); rot(b); break;
        case RRA: rrot(a); break; case RRB: rrot(b); break; case RRR: rrot(a); rrot(b); break;
        default: break;
        }
    }
    bool sorted() const {
        if (!b.empty()) return false;
        for (size_t i = 1; i < a.size(); i++) if (a[i - 1] > a[i]) return false;
        return true;
    }
};

static inline bool verify_ops(const std::vector<int> &init, const std::vector<Op> &ops) {
    Sim s; for (int x : init) s.a.push_back(x);
    for (Op o : ops) s.apply(o);
    return s.sorted();
}

// Parse argv into values; returns false on error (non-int, overflow, duplicate).
static inline bool parse_args(int argc, char **argv, int first, std::vector<long long> &vals) {
    for (int i = first; i < argc; i++) {
        const char *s = argv[i];
        while (*s) {
            while (*s == ' ' || *s == '\t') s++;
            if (!*s) break;
            const char *st = s; char *end;
            if (*s == '+' || *s == '-') s++;
            if (!(*s >= '0' && *s <= '9')) return false;
            errno = 0; long long v = strtoll(st, &end, 10);
            if (errno || v < INT_MIN || v > INT_MAX) return false;
            if (*end && *end != ' ' && *end != '\t') return false;
            vals.push_back(v); s = end;
        }
    }
    std::vector<long long> c = vals; std::sort(c.begin(), c.end());
    for (size_t i = 1; i < c.size(); i++) if (c[i] == c[i - 1]) return false;
    return true;
}

static inline std::vector<int> to_ranks(const std::vector<long long> &vals) {
    std::vector<long long> c = vals; std::sort(c.begin(), c.end());
    std::vector<int> r(vals.size());
    for (size_t i = 0; i < vals.size(); i++) r[i] = (int)(std::lower_bound(c.begin(), c.end(), vals[i]) - c.begin());
    return r;
}

static inline void print_ops(const std::vector<Op> &ops) {
    std::string out; out.reserve(ops.size() * 4);
    for (Op o : ops) { out += OPNAME[o]; out += '\n'; }
    fwrite(out.data(), 1, out.size(), stdout);
}

struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed = 88172645463325252ull) : s(seed ? seed : 88172645463325252ull) {}
    inline uint64_t next() { s ^= s << 7; s ^= s >> 9; return s; }
    inline int below(int n) { return (int)((next() >> 11) % (uint64_t)n); }
    inline double uni() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
};

static inline int comb_cost_simple(int fa, int sa, int fb, int sb) {
    int ra = fa, rra = fa ? sa - fa : 0, rb = fb, rrb = fb ? sb - fb : 0;
    return std::min(std::min(std::max(ra, rb), std::max(rra, rrb)), std::min(ra + rrb, rra + rb));
}
