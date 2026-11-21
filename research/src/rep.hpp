// Solution representation (query list + AOrder/BOrder keys), decoder, converter.
//
// Representation:
//   inK[v]  : v never moves (K must be a circularly increasing subsequence of the input).
//   q       : sequence of queries A2B(v)/B2A(v) for every v not in K, A2B(v) before B2A(v).
//   ik[v]   : AOrder key of v in its initial state (position in initial A).
//   tk[v]   : AOrder key of v after returning to A (target). For K values tk == ik.
//   bk[v]   : BOrder key of v.
// Stacks are kept cyclically sorted by key; an insertion goes to the unique place that keeps
// the cyclic order. Targets are cyclically sorted by value => final A is sorted after rotation.
#pragma once
#include "common.hpp"

typedef long long ll;
struct Query { int16_t v; uint8_t t; };  // t: 0 = A2B, 1 = B2A

static const ll KEY_BASE = 1LL << 40;
static const ll KEY_GAP = 1LL << 22;

struct Rep {
    int n = 0;
    std::vector<int> a0;
    std::vector<char> inK;
    std::vector<Query> q;
    std::vector<ll> ik, tk, bk;
    std::vector<char> present;  // values currently part of the instance (IG destroy removes values)
    int minv = 0;               // smallest present value (final rotation target)
    void fix_minv() { minv = 0; while (minv < n && !present[minv]) minv++; }
};

static inline int comb_cost(int fa, int sa, int fb, int sb, int *mode = nullptr) {
    int ra = fa, rra = fa ? sa - fa : 0, rb = fb, rrb = fb ? sb - fb : 0;
    int c0 = std::max(ra, rb), c1 = std::max(rra, rrb), c2 = ra + rrb, c3 = rra + rb;
    int best = c0, m = 0;
    if (c1 < best) best = c1, m = 1;
    if (c2 < best) best = c2, m = 2;
    if (c3 < best) best = c3, m = 3;
    if (mode) *mode = m;
    return best;
}

// Stack state: sorted encoded entries (key*1024+v) and index of top.
struct St {
    std::vector<ll> A, B;
    int ta = 0, tb = 0;
};

static inline int topv(const std::vector<ll> &s, int t) { return s.empty() ? -1 : (int)(s[t] & 1023); }

// One decoding step. Returns cost (rotations + push). Optionally reports rotation info.
struct StepInfo { int fa, sa, fb, sb, mode; };

static inline int do_step(St &s, const Rep &r, Query qu, StepInfo *info = nullptr) {
    int v = qu.v;
    std::vector<ll> &S = qu.t == 0 ? s.A : s.B;  // source
    std::vector<ll> &D = qu.t == 0 ? s.B : s.A;  // destination
    int &ts = qu.t == 0 ? s.ta : s.tb;
    int &td = qu.t == 0 ? s.tb : s.ta;
    ll es = (qu.t == 0 ? r.ik[v] : r.bk[v]) * 1024 + v;
    ll ed = (qu.t == 0 ? r.bk[v] : r.tk[v]) * 1024 + v;
    int ss = S.size(), sd = D.size();
    int i = (int)(std::lower_bound(S.begin(), S.end(), es) - S.begin());
    int fs = ((i - ts) % ss + ss) % ss;
    int j = (int)(std::lower_bound(D.begin(), D.end(), ed) - D.begin());
    int fd = 0;
    if (sd > 0) { int jj = j == sd ? 0 : j; fd = ((jj - td) % sd + sd) % sd; }
    int fa = qu.t == 0 ? fs : fd, sa = qu.t == 0 ? ss : sd, fb = qu.t == 0 ? fd : fs, sb = qu.t == 0 ? sd : ss;
    int mode;
    int c = comb_cost(fa, sa, fb, sb, &mode) + 1;
    if (info) *info = {fa, sa, fb, sb, mode};
    S.erase(S.begin() + i);
    ts = S.empty() ? 0 : (i == (int)S.size() ? 0 : i);
    D.insert(D.begin() + j, ed);
    td = j;
    return c;
}

static inline int final_rot(const St &s, const Rep &r, int *dir = nullptr) {
    int sa = s.A.size(); if (sa <= 1 || r.minv >= r.n) { if (dir) *dir = 0; return 0; }
    ll e = r.tk[r.minv] * 1024 + r.minv;
    int i = (int)(std::lower_bound(s.A.begin(), s.A.end(), e) - s.A.begin());
    int f = ((i - s.ta) % sa + sa) % sa;
    if (dir) *dir = f <= sa - f ? 0 : 1;
    return std::min(f, sa - f);
}

static inline St initial_state(const Rep &r) {
    St s;
    for (int x : r.a0) if (r.present[x]) s.A.push_back(r.ik[x] * 1024 + x);
    // a0 order is AOrder order (keys increasing along a0 from a0[0]) only up to rotation; sort and find top.
    std::sort(s.A.begin(), s.A.end());
    int first = -1; for (int x : r.a0) if (r.present[x]) { first = x; break; }
    s.ta = 0;
    if (first >= 0) {
        ll e0 = r.ik[first] * 1024 + first;
        s.ta = (int)(std::lower_bound(s.A.begin(), s.A.end(), e0) - s.A.begin());
    }
    s.tb = 0;
    return s;
}

static inline int full_cost(const Rep &r) {
    St s = initial_state(r);
    int c = 0;
    for (auto qu : r.q) c += do_step(s, r, qu);
    return c + final_rot(s, r);
}

static inline std::vector<Op> emit(const Rep &r) {
    std::vector<Op> ops;
    St s = initial_state(r);
    for (auto qu : r.q) {
        StepInfo in;
        do_step(s, r, qu, &in);
        int ra = in.fa, rra = in.fa ? in.sa - in.fa : 0, rb = in.fb, rrb = in.fb ? in.sb - in.fb : 0;
        auto put = [&](Op o, int k) { for (int i = 0; i < k; i++) ops.push_back(o); };
        switch (in.mode) {
        case 0: put(RR, std::min(ra, rb)); put(RA, ra - std::min(ra, rb)); put(RB, rb - std::min(ra, rb)); break;
        case 1: put(RRR, std::min(rra, rrb)); put(RRA, rra - std::min(rra, rrb)); put(RRB, rrb - std::min(rra, rrb)); break;
        case 2: put(RA, ra); put(RRB, rrb); break;
        case 3: put(RRA, rra); put(RB, rb); break;
        }
        ops.push_back(qu.t == 0 ? PB : PA);
    }
    int dir, c = final_rot(s, r, &dir);
    for (int i = 0; i < c; i++) ops.push_back(dir == 0 ? RA : RRA);
    return ops;
}

// Convert a swap-free command sequence (each value pushed to B at most once) into a Rep.
// Uses ghost lists: removed entries stay in place, so a single global cyclic order exists.
static inline bool convert(const std::vector<int> &a0, const std::vector<Op> &ops, Rep &r) {
    int n = a0.size();
    struct E { int v; bool live; bool target; };
    std::vector<E> A, B;
    for (int x : a0) A.push_back({x, true, false});
    int ta = 0, tb = 0, liveB = 0;
    std::vector<int> trips(n, 0);
    r = Rep(); r.n = n; r.a0 = a0; r.inK.assign(n, 1); r.present.assign(n, 1); r.minv = 0; r.ik.assign(n, 0); r.tk.assign(n, 0); r.bk.assign(n, 0);
    auto nextlive = [](std::vector<E> &L, int t, int dir) {
        int m = L.size();
        for (int k = 1; k <= m; k++) { int i = ((t + dir * k) % m + m) % m; if (L[i].live) return i; }
        return t;
    };
    int liveA = n;
    for (Op o : ops) {
        switch (o) {
        case RA: case RR: case RRA: case RRR: case RB: case RRB: break;
        case PB: case PA: break;
        default: return false;
        }
        if (o == RA || o == RR) { if (liveA > 1) ta = nextlive(A, ta, 1); }
        if (o == RRA || o == RRR) { if (liveA > 1) ta = nextlive(A, ta, -1); }
        if (o == RB || o == RR) { if (liveB > 1) tb = nextlive(B, tb, 1); }
        if (o == RRB || o == RRR) { if (liveB > 1) tb = nextlive(B, tb, -1); }
        if (o == PB) {
            if (!liveA) return false;
            int v = A[ta].v; if (A[ta].target || trips[v]) return false;
            trips[v]++; r.inK[v] = 0; r.q.push_back({(int16_t)v, 0});
            A[ta].live = false; liveA--;
            if (liveA) ta = nextlive(A, ta, 1);
            int pos = liveB ? tb : (int)B.size();
            B.insert(B.begin() + pos, {v, true, false}); tb = pos; liveB++;
        }
        if (o == PA) {
            if (!liveB) return false;
            int v = B[tb].v; B[tb].live = false; liveB--;
            r.q.push_back({(int16_t)v, 1});
            if (liveB) tb = nextlive(B, tb, 1);
            int pos = ta;  // insert before current top (A never empty here if K nonempty)
            if (!liveA) pos = A.size();
            A.insert(A.begin() + pos, {v, true, true}); ta = pos; liveA++;
        }
    }
    if (liveB) return false;
    for (int i = 0; i < (int)A.size(); i++) {
        ll key = KEY_BASE + (ll)i * KEY_GAP;
        if (A[i].target) r.tk[A[i].v] = key; else r.ik[A[i].v] = key;
    }
    for (int v = 0; v < n; v++) if (r.inK[v]) r.tk[v] = r.ik[v];
    for (int i = 0; i < (int)B.size(); i++) r.bk[B[i].v] = KEY_BASE + (ll)i * KEY_GAP;
    return true;
}
