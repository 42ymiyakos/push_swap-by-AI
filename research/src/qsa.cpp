// QSA: query-representation simulated annealing solver for push_swap.
// Build: g++ -O2 -std=c++17 -o qsa qsa.cpp
#include "common.hpp"
#include "rep.hpp"
#include "construct.hpp"
#include "exact.hpp"
#include <chrono>
#include <cmath>

static double now_sec() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

struct Params {
    double time = 5.0;
    uint64_t seed = 1;
    int chunks = -1;
    double T0 = 3.0, T1 = 0.2;
    double p_shift = 0.6, p_border = 0.25, p_aorder = 0.1, p_kadd = 0.025, p_krem = 0.025;
    int shift_w = 12, border_w = 6;
    int verbose = 0;
    int check = 0;
};

struct AE { ll key; int v; int kind; };  // kind 0 init, 1 target, 2 both (K)

struct Opt {
    Rep r;
    Params P;
    Rng rng;
    static const int C = 32;
    std::vector<St> cps;
    std::vector<int> cum;
    std::vector<int16_t> rtA, rtB;
    int total = 0;
    std::vector<int> pos[2];
    std::vector<int> bl;   // non-K values sorted by bk
    std::vector<AE> al;    // all A entries sorted by key
    St scr;
    long long n_eval = 0, n_acc = 0, steps = 0;
    long long cnt_try[5] = {0}, cnt_acc[5] = {0};

    void recompute_pos() {
        pos[0].assign(r.n, -1); pos[1].assign(r.n, -1);
        for (int i = 0; i < (int)r.q.size(); i++) pos[r.q[i].t][r.q[i].v] = i;
    }
    void rebuild_lists() {
        bl.clear();
        for (int v = 0; v < r.n; v++) if (!r.inK[v]) bl.push_back(v);
        std::sort(bl.begin(), bl.end(), [&](int x, int y) { return r.bk[x] < r.bk[y]; });
        al.clear();
        for (int v = 0; v < r.n; v++) {
            if (r.inK[v]) al.push_back({r.ik[v], v, 2});
            else { al.push_back({r.ik[v], v, 0}); al.push_back({r.tk[v], v, 1}); }
        }
        std::sort(al.begin(), al.end(), [](const AE &x, const AE &y) { return x.key < y.key; });
    }
    void rebuild(int i0) {
        int Q = r.q.size();
        int c = i0 / C;
        if (c >= (int)cps.size()) c = (int)cps.size() - 1;
        if (c < 0 || cps.empty()) { c = 0; cps.assign(1, initial_state(r)); }
        if (c == 0) cps[0] = initial_state(r);
        cps.resize(Q / C + 1);
        cum.resize(Q + 1); rtA.resize(Q + 1); rtB.resize(Q + 1);
        St s = cps[c];
        if (c == 0) cum[0] = 0;
        for (int i = c * C; i < Q; i++) {
            if (i % C == 0) cps[i / C] = s;
            rtA[i] = topv(s.A, s.ta); rtB[i] = topv(s.B, s.tb);
            cum[i + 1] = cum[i] + do_step(s, r, r.q[i]);
        }
        rtA[Q] = topv(s.A, s.ta); rtB[Q] = topv(s.B, s.tb);
        total = cum[Q] + final_rot(s, r);
    }
    int evaluate(int i0, int stopFrom, int off) {
        n_eval++;
        int Q = r.q.size();
        int c = i0 / C;
        scr = cps[c];
        int cost = cum[c * C];
        for (int i = c * C; i < Q; i++) {
            if (i >= stopFrom) {
                int oi = i - off;
                if (topv(scr.A, scr.ta) == rtA[oi] && topv(scr.B, scr.tb) == rtB[oi]) return cost + total - cum[oi];
            }
            cost += do_step(scr, r, r.q[i]);
        }
        return cost + final_rot(scr, r);
    }
    bool accept(int nc, double T) {
        int d = nc - total;
        if (d <= 0) return true;
        return rng.uni() < std::exp(-d / T);
    }

    // ---- move: shift a query ----
    void do_shift(int i, int p) {
        if (i < p) std::rotate(r.q.begin() + i, r.q.begin() + i + 1, r.q.begin() + p + 1);
        else std::rotate(r.q.begin() + p, r.q.begin() + i, r.q.begin() + i + 1);
        int lo = std::min(i, p), hi = std::max(i, p);
        for (int k = lo; k <= hi; k++) pos[r.q[k].t][r.q[k].v] = k;
    }
    bool mv_shift(double T) {
        int Q = r.q.size(); if (Q < 2) return false;
        int i = rng.below(Q);
        int w = P.shift_w;
        int d = 1 + rng.below(w);
        if (rng.below(8) == 0) d = 1 + rng.below(std::max(1, Q / 4));
        if (rng.below(2)) d = -d;
        int p = std::clamp(i + d, 0, Q - 1);
        if (p == i) return false;
        Query x = r.q[i];
        if (x.t == 0) { int j = pos[1][x.v]; if (p >= j) return false; }
        else { int j = pos[0][x.v]; if (p <= j) return false; }
        do_shift(i, p);
        int nc = evaluate(std::min(i, p), std::max(i, p) + 1, 0);
        if (accept(nc, T)) { rebuild(std::min(i, p)); return true; }
        do_shift(p, i);
        return false;
    }

    // ---- move: BOrder ----
    void respace_b() {
        for (int i = 0; i < (int)bl.size(); i++) r.bk[bl[i]] = KEY_BASE + (ll)i * KEY_GAP;
        rebuild(0);
    }
    bool mv_border(double T) {
        int m = bl.size(); if (m < 3) return false;
        int r0 = rng.below(m), v = bl[r0];
        int L = m - 1;  // others
        int k0 = (r0 - 1 + L) % L;  // slot after others[k0]  (others index == bl index skipping v)
        int d = 1 + rng.below(P.border_w);
        if (rng.below(8) == 0) d = 1 + rng.below(L);
        if (rng.below(2)) d = -d;
        int k = ((k0 + d) % L + L) % L;
        if (k == k0) return false;
        auto oth = [&](int idx) { return bl[idx < r0 ? idx : idx + 1]; };
        ll nk;
        if (k < L - 1) {
            ll a = r.bk[oth(k)], b = r.bk[oth(k + 1)];
            if (b - a < 2) { respace_b(); return false; }
            nk = a + (b - a) / 2;
        } else nk = r.bk[oth(L - 1)] + KEY_GAP;
        ll ok = r.bk[v];
        r.bk[v] = nk;
        int nc = evaluate(pos[0][v], pos[1][v] + 1, 0);
        if (accept(nc, T)) {
            bl.erase(bl.begin() + r0);
            bl.insert(std::upper_bound(bl.begin(), bl.end(), v, [&](int x, int y) { return r.bk[x] < r.bk[y]; }), v);
            rebuild(pos[0][v]);
            return true;
        }
        r.bk[v] = ok;
        return false;
    }

    // ---- move: AOrder target relocation ----
    int al_index(ll key) {
        return (int)(std::lower_bound(al.begin(), al.end(), key, [](const AE &e, ll k) { return e.key < k; }) - al.begin());
    }
    bool mv_aorder(double T) {
        int m = bl.size(); if (m < 1 || r.n < 3) return false;
        int v = bl[rng.below(m)];
        int pv = (v - 1 + r.n) % r.n, sv = (v + 1) % r.n;
        int L = al.size();
        int iv = al_index(r.tk[v]);
        int ip = al_index(r.tk[pv]), is = al_index(r.tk[sv]);
        // indices in list without v's target entry
        auto adj = [&](int x) { return x > iv ? x - 1 : x; };
        ip = adj(ip); is = adj(is);
        int L2 = L - 1;
        int cnt = ((is - ip) % L2 + L2) % L2;
        if (cnt <= 1) return false;
        int kcur = ((iv - 1) % L2 + L2) % L2;  // slot after entry kcur (in reduced list)
        int k = (ip + rng.below(cnt)) % L2;
        if (k == kcur) return false;
        auto ent = [&](int idx) { return al[idx < iv ? idx : idx + 1]; };
        ll nk;
        if (k < L2 - 1) {
            ll a = ent(k).key, b = ent(k + 1).key;
            if (b - a < 2) { respace_a(); return false; }
            nk = a + (b - a) / 2;
        } else nk = ent(L2 - 1).key + KEY_GAP;
        ll ok = r.tk[v];
        r.tk[v] = nk;
        int nc = evaluate(pos[1][v], INT_MAX, 0);
        if (accept(nc, T)) {
            AE e = al[iv]; al.erase(al.begin() + iv); e.key = nk;
            al.insert(std::lower_bound(al.begin(), al.end(), nk, [](const AE &x, ll kk) { return x.key < kk; }), e);
            rebuild(pos[1][v]);
            return true;
        }
        r.tk[v] = ok;
        return false;
    }
    void respace_a() {
        for (int i = 0; i < (int)al.size(); i++) {
            ll key = KEY_BASE + (ll)i * KEY_GAP; al[i].key = key;
            if (al[i].kind == 0) r.ik[al[i].v] = key;
            else if (al[i].kind == 1) r.tk[al[i].v] = key;
            else r.ik[al[i].v] = r.tk[al[i].v] = key;
        }
        rebuild(0);
    }

    static bool between(ll a, ll x, ll b) { return a < b ? (a < x && x < b) : (x > a || x < b); }

    // ---- move: add v to K (stop moving it) ----
    bool mv_kadd(double T) {
        int m = bl.size(); if (m < 1 || r.n < 3) return false;
        int v = bl[rng.below(m)];
        int pv = (v - 1 + r.n) % r.n, sv = (v + 1) % r.n;
        if (!between(r.tk[pv], r.ik[v], r.tk[sv])) return false;
        Rep save = r;
        int i0 = pos[0][v], i1 = pos[1][v];
        r.q.erase(r.q.begin() + i1); r.q.erase(r.q.begin() + i0);
        r.inK[v] = 1; r.tk[v] = r.ik[v];
        int nc = evaluate(i0, INT_MAX, 0);
        if (accept(nc, T)) { recompute_pos(); rebuild_lists(); rebuild(i0); return true; }
        r = save;
        return false;
    }
    // ---- move: remove v from K (send it through B) ----
    bool mv_krem(double T) {
        if (r.n < 3) return false;
        int kc = 0; for (int v = 0; v < r.n; v++) kc += r.inK[v];
        if (kc <= 1) return false;
        int v; do v = rng.below(r.n); while (!r.inK[v]);
        int Q = r.q.size();
        // A2B position: near A2B of an initial-A neighbour (cheap rotation), else random
        int idx0 = 0; for (int i = 0; i < r.n; i++) if (r.a0[i] == v) idx0 = i;
        int p = rng.below(Q + 1);
        for (int t = 0; t < 4; t++) {
            int nb = r.a0[(idx0 + (rng.below(2) ? 1 : r.n - 1) * (1 + rng.below(2))) % r.n];
            if (!r.inK[nb]) { p = pos[0][nb] + rng.below(2); break; }
        }
        // B2A position: next to B2A of a value neighbour
        int p2 = -1; int nbv = -1;
        for (int t = 0; t < 6; t++) {
            int w = (v + (rng.below(2) ? 1 : -1) * (1 + rng.below(3)) + r.n) % r.n;
            if (!r.inK[w] && pos[1][w] + 1 > p) { p2 = pos[1][w] + rng.below(2); nbv = w; break; }
        }
        if (p2 < 0) p2 = p + 1 + rng.below(std::max(1, Q - p));
        if (p2 <= p) p2 = p + 1;
        Rep save = r;
        r.q.insert(r.q.begin() + p, {(int16_t)v, 0});
        p2 = std::min(p2 + 1, (int)r.q.size());
        r.q.insert(r.q.begin() + p2, {(int16_t)v, 1});
        r.inK[v] = 0;
        // bk next to the value-neighbour (or random)
        if (nbv >= 0) {
            r.bk[v] = r.bk[nbv] + (rng.below(2) ? 1 : -1) * (KEY_GAP / 64 + rng.below(1000));
        } else if (!bl.empty()) r.bk[v] = r.bk[bl[rng.below(bl.size())]] + 1 + rng.below(1000);
        // tk right after ik
        int ia = al_index(r.ik[v]);
        ll nxt = ia + 1 < (int)al.size() ? al[ia + 1].key : r.ik[v] + 2 * KEY_GAP;
        if (nxt - r.ik[v] < 2) { r = save; respace_a(); return false; }
        r.tk[v] = r.ik[v] + (nxt - r.ik[v]) / 2;
        // uniqueness of bk
        for (int u : bl) if (r.bk[u] == r.bk[v]) { r = save; return false; }
        int nc = evaluate(p, INT_MAX, 0);
        if (accept(nc, T)) { recompute_pos(); rebuild_lists(); rebuild(p); return true; }
        r = save;
        return false;
    }

    void init(const Rep &rep) {
        r = rep; recompute_pos(); rebuild_lists(); cps.clear(); rebuild(0);
    }

    Rep best; int best_cost = INT_MAX;

    void run(double tlimit) {
        double t0 = now_sec();
        best = r; best_cost = total;
        double T = P.T0, prog = 0;
        double pc[5] = {P.p_shift, P.p_border, P.p_aorder, P.p_kadd, P.p_krem};
        double ps = 0; for (double x : pc) ps += x;
        for (long long it = 0;; it++) {
            if ((it & 255) == 0) {
                prog = (now_sec() - t0) / tlimit;
                if (prog >= 1) break;
                T = P.T0 * std::pow(P.T1 / P.T0, prog);
            }
            double u = rng.uni() * ps; int mt = 0;
            while (mt < 4 && u >= pc[mt]) { u -= pc[mt]; mt++; }
            bool ok = false;
            cnt_try[mt]++;
            switch (mt) {
            case 0: ok = mv_shift(T); break;
            case 1: ok = mv_border(T); break;
            case 2: ok = mv_aorder(T); break;
            case 3: ok = mv_kadd(T); break;
            case 4: ok = mv_krem(T); break;
            }
            if (ok) {
                cnt_acc[mt]++;
                n_acc++;
                if (P.check && full_cost(r) != total) { fprintf(stderr, "MISMATCH mt=%d %d %d\n", mt, full_cost(r), total); exit(3); }
                if (total < best_cost) { best_cost = total; best = r; }
            }
            steps++;
        }
    }
};

int main(int argc, char **argv) {
    Params P;
    int first = 1;
    for (; first < argc && argv[first][0] == '-' && argv[first][1] == '-'; first++) {
        std::string a = argv[first];
        auto val = [&](const char *k) -> const char * { size_t L = strlen(k); return a.compare(0, L, k) == 0 ? a.c_str() + L : nullptr; };
        if (auto x = val("--time=")) P.time = atof(x);
        else if (auto x = val("--seed=")) P.seed = strtoull(x, 0, 10);
        else if (auto x = val("--chunks=")) P.chunks = atoi(x);
        else if (auto x = val("--T0=")) P.T0 = atof(x);
        else if (auto x = val("--T1=")) P.T1 = atof(x);
        else if (auto x = val("--pshift=")) P.p_shift = atof(x);
        else if (auto x = val("--pborder=")) P.p_border = atof(x);
        else if (auto x = val("--paorder=")) P.p_aorder = atof(x);
        else if (auto x = val("--pkadd=")) P.p_kadd = atof(x);
        else if (auto x = val("--pkrem=")) P.p_krem = atof(x);
        else if (auto x = val("--shiftw=")) P.shift_w = atoi(x);
        else if (auto x = val("--borderw=")) P.border_w = atoi(x);
        else if (auto x = val("--verbose=")) P.verbose = atoi(x);
        else if (auto x = val("--check=")) P.check = atoi(x);
        else { fprintf(stderr, "Error\n"); return 1; }
    }
    std::vector<long long> vals;
    if (!parse_args(argc, argv, first, vals)) { fprintf(stderr, "Error\n"); return 1; }
    int n = vals.size();
    if (n == 0) return 0;
    std::vector<int> a0 = to_ranks(vals);
    bool sorted = true; for (int i = 0; i < n; i++) if (a0[i] != i) sorted = false;
    if (sorted) return 0;
    if (n <= 7) { auto ops = exact::solve(a0); print_ops(ops); return 0; }

    Rng rng(P.seed * 0x9E3779B97F4A7C15ull + 12345);
    int m = P.chunks > 0 ? P.chunks : (n <= 150 ? 2 : (n <= 300 ? 4 : 6));
    int lis_len;
    std::vector<char> K = circular_lis(a0, rng, &lis_len);
    std::vector<Op> base = chunk_construct(a0, K, m);
    if (!verify_ops(a0, base)) { fprintf(stderr, "internal: constructor invalid\n"); return 2; }
    Rep rep;
    if (!convert(a0, base, rep)) { fprintf(stderr, "internal: convert failed\n"); print_ops(base); return 0; }
    int c0 = full_cost(rep);
    Opt opt; opt.P = P; opt.rng = Rng(P.seed + 777);
    opt.init(rep);
    double t0 = now_sec();
    opt.run(P.time);
    std::vector<Op> ops = emit(opt.best);
    if ((int)ops.size() != opt.best_cost || !verify_ops(a0, ops)) {
        fprintf(stderr, "internal: emit mismatch %zu vs %d\n", ops.size(), opt.best_cost);
        ops = base;
    }
    if (P.verbose) {
        fprintf(stderr, "n=%d lis=%d constructor=%zu rep=%d final=%zu time=%.2f steps=%lld acc=%lld\n", n, lis_len,
                base.size(), c0, ops.size(), now_sec() - t0, opt.steps, opt.n_acc);
        const char *nm[5] = {"shift", "border", "aorder", "kadd", "krem"};
        for (int i = 0; i < 5; i++) fprintf(stderr, "  %s try=%lld acc=%lld\n", nm[i], opt.cnt_try[i], opt.cnt_acc[i]);
    }
    print_ops(ops);
    return 0;
}
