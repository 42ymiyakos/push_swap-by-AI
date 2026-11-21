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
struct Query { int16_t v; uint8_t t; };  // t: 0 = A2B, 1 = B2A, 2 = A2A (move inside A with sa)

static const ll KEY_BASE = 1LL << 40;
static const ll KEY_GAP = 1LL << 22;

struct Rep {
    int n = 0;
    std::vector<int> a0;
    std::vector<char> inK;
    std::vector<Query> q;
    std::vector<ll> ik, tk, bk;
    std::vector<char> isA2A;    // v is handled by an A2A query (never visits B)
    std::vector<char> ph;       // placeholder: query-destroyed value sitting in A at its init key
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

// A2A: element u at index src of A (size sa, top index ta) moves to the gap given by dstRaw
// (lower_bound of its target key in A, u included). Executed with sa and rotations:
//   forward d steps : rotate u to top, then sa,(ra,sa)^(d-1)      -> new top = last jumped element
//   backward d steps: rotate pred(u) to top, then sa,(rra,sa)^(d-1) -> new top = u
struct A2AInfo { int cost, fwd, steps, rot, rotfwd, topWo; };  // topWo: index (in A without u) of new top if fwd
static inline A2AInfo a2a_cost(int sa, int ta, int src, int dstRaw) {
    A2AInfo r{0, 1, 0, 0, 1, -1};
    int m = sa - 1;
    if (m <= 0) return r;
    int rt = dstRaw - (dstRaw > src ? 1 : 0);
    int df = ((rt - src) % m + m) % m, db = (m - df) % m;
    if (df == 0) { r.topWo = -2; return r; }  // no move; top unchanged (caller keeps element)
    int f1 = ((src - ta) % sa + sa) % sa, r1 = std::min(f1, sa - f1);
    int pi = (src - 1 + sa) % sa;
    int f2 = ((pi - ta) % sa + sa) % sa, r2 = std::min(f2, sa - f2);
    int cf = r1 + 2 * df - 1, cb = r2 + 2 * db - 1;
    if (cf <= cb) { r = {cf, 1, df, r1, f1 <= sa - f1, (src + df - 1) % m}; }
    else { r = {cb, 0, db, r2, f2 <= sa - f2, -1}; }
    return r;
}

// Push cost with pending A rotations (pr = ra, prr = rra) from a preceding A2A run that can be
// paired with this push's B rotations (rr / rrr). mode: combo; *k: number of paired rotations.
static inline int comb_sync(int fa, int sa, int fb, int sb, int pr, int prr, int *mode = nullptr, int *kk = nullptr) {
    int ra = fa, rra = fa ? sa - fa : 0, rb = fb, rrb = fb ? sb - fb : 0;
    int k0 = std::min(pr, std::max(0, rb - ra)), k1 = std::min(prr, std::max(0, rrb - rra));
    int k2 = std::min(prr, rrb), k3 = std::min(pr, rb);
    int c0 = std::max(ra, rb) - k0, c1 = std::max(rra, rrb) - k1, c2 = ra + rrb - k2, c3 = rra + rb - k3;
    int best = c0, m = 0, k = k0;
    if (c1 < best) best = c1, m = 1, k = k1;
    if (c2 < best) best = c2, m = 2, k = k2;
    if (c3 < best) best = c3, m = 3, k = k3;
    if (mode) *mode = m;
    if (kk) *kk = k;
    return best;
}

// Stack state: sorted encoded entries (key*1024+v) and index of top.
struct St {
    std::vector<ll> A, B;
    int ta = 0, tb = 0;
    int pr = 0, prr = 0;  // pending A rotations of the current A2A run (sync with next push)
};

static inline int topv(const std::vector<ll> &s, int t) { return s.empty() ? -1 : (int)(s[t] & 1023); }

// One decoding step. Returns cost (rotations + push). Optionally reports rotation info.
struct StepInfo { int fa, sa, fb, sb, mode, k; };

static inline int do_step_a2a(St &s, const Rep &r, Query qu, A2AInfo *info = nullptr) {
    int v = qu.v;
    std::vector<ll> &A = s.A;
    ll es = r.ik[v] * 1024 + v, ed = r.tk[v] * 1024 + v;
    int sa = A.size();
    int src = (int)(std::lower_bound(A.begin(), A.end(), es) - A.begin());
    int dr = (int)(std::lower_bound(A.begin(), A.end(), ed) - A.begin());
    A2AInfo in = a2a_cost(sa, s.ta, src, dr);
    if (info) *info = in;
    if (in.topWo != -2 && sa > 1) {
        if (in.rotfwd) s.pr += in.rot; else s.prr += in.rot;
        if (in.fwd) s.pr += in.steps - 1; else s.prr += in.steps - 1;
    }
    ll topKey = A.empty() ? 0 : A[s.ta];
    A.erase(A.begin() + src);
    int j = (int)(std::lower_bound(A.begin(), A.end(), ed) - A.begin());
    ll newTop;
    if (in.topWo == -2 || in.topWo == -1 && sa <= 1) newTop = topKey == es ? ed : topKey;
    else if (in.fwd) newTop = A[in.topWo];
    else newTop = ed;
    A.insert(A.begin() + j, ed);
    s.ta = (int)(std::lower_bound(A.begin(), A.end(), newTop) - A.begin());
    return in.cost;
}

static inline int do_step(St &s, const Rep &r, Query qu, StepInfo *info = nullptr) {
    if (qu.t == 2) { if (info) *info = {0, 0, 0, 0, -1, 0}; return do_step_a2a(s, r, qu); }
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
    int mode, kk;
    int c = comb_sync(fa, sa, fb, sb, s.pr, s.prr, &mode, &kk) + 1;
    if (info) *info = {fa, sa, fb, sb, mode, kk};
    s.pr = s.prr = 0;
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
    size_t runStart = 0;  // start index of the current A2A run in ops
    for (auto qu : r.q) {
        if (qu.t == 2) {
            A2AInfo a; int sa = s.A.size();
            do_step_a2a(s, r, qu, &a);
            if (a.topWo == -2 || sa <= 1) continue;
            for (int k = 0; k < a.rot; k++) ops.push_back(a.rotfwd ? RA : RRA);
            ops.push_back(SA);
            for (int k = 1; k < a.steps; k++) { ops.push_back(a.fwd ? RA : RRA); ops.push_back(SA); }
            continue;
        }
        StepInfo in;
        do_step(s, r, qu, &in);
        // pair k pending A rotations of the run with B rotations of this push
        if (in.k > 0) {
            Op want = (in.mode == 0 || in.mode == 3) ? RA : RRA;  // B dir fwd for modes 0,3
            Op merged = want == RA ? RR : RRR;
            int k = in.k;
            for (size_t x = runStart; x < ops.size() && k > 0; x++) if (ops[x] == want) { ops[x] = merged; k--; }
        }
        int kB = in.k;
        int ra = in.fa, rra = in.fa ? in.sa - in.fa : 0, rb = in.fb, rrb = in.fb ? in.sb - in.fb : 0;
        auto put = [&](Op o, int k) { for (int i = 0; i < k; i++) ops.push_back(o); };
        switch (in.mode) {
        case 0: put(RR, std::min(ra, rb)); put(RA, ra - std::min(ra, rb)); put(RB, rb - std::min(ra, rb) - kB); break;
        case 1: put(RRR, std::min(rra, rrb)); put(RRA, rra - std::min(rra, rrb)); put(RRB, rrb - std::min(rra, rrb) - kB); break;
        case 2: put(RA, ra); put(RRB, rrb - kB); break;
        case 3: put(RRA, rra); put(RB, rb - kB); break;
        }
        ops.push_back(qu.t == 0 ? PB : PA);
        runStart = ops.size();
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
    r = Rep(); r.n = n; r.a0 = a0; r.inK.assign(n, 1); r.isA2A.assign(n, 0); r.ph.assign(n, 0); r.present.assign(n, 1); r.minv = 0; r.ik.assign(n, 0); r.tk.assign(n, 0); r.bk.assign(n, 0);
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
        case SA: {  // swap of two init-state values in A -> A2A query (top value jumps one step forward)
            if (liveA < 2) continue;
            int x = A[ta].v; if (A[ta].target || trips[x]) return false;
            int ny = nextlive(A, ta, 1);
            trips[x]++; r.inK[x] = 0; r.isA2A[x] = 1; r.q.push_back({(int16_t)x, 2});
            A[ta].live = false;
            A.insert(A.begin() + ny + 1, {x, true, true});
            if (ny + 1 <= ta) ta++;
            ta = ny < ta ? ny : ny;  // y is the new top
            continue;
        }
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
