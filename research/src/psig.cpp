// PSIG: push_swap Iterated Greedy with exact DP re-insertion + LAHC acceptance.
// Build: g++ -O2 -std=c++17 -o psig psig.cpp
#include "common.hpp"
#include "rep.hpp"
#include "construct.hpp"
#include "exact.hpp"
#include "ig.hpp"
#include "peephole.hpp"
#include <chrono>
#include <cmath>

static double now_sec() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

struct Params {
    double time = 5.0;
    long long iters = -1;
    uint64_t seed = 1;
    int chunks = -1;
    int bwin = 4;        // value window for B slot candidates
    int ball = 1;        // 1: all B slots are candidates
    int tmode = 1, t0 = 0, topb = 3;
    int twophase = 0;
    double pqd = 0;
    std::vector<int> peepK{};
    int restarts = 1; double rfrac = 0.5; int kinit = 1;
    int brand = 2, ginit = 0, sweep = 0, a2a = 0, a2amax = 3;
    int bw = 1000, bws = 20, bc = 2, nlis = 50, ktop = 3;
    int seglen = 8;      // max query-segment length for destroy
    int rmax = 6;        // max random destroy count
    int vwin = 6;        // max value-range destroy width
    int lahc = 1;      // LAHC history length
    int verbose = 0;
    int check = 0;
};

struct IG {
    Rep cur; int curCost = 0;
    Rep best; int bestCost = INT_MAX;
    Params P; Rng rng;
    Reinserter R;
    std::vector<int> vcost;
    long long nMis = 0;
    Reinserter::Best lastInfo; int lastQ = 0; bool lastVBF = false;
    long long iters = 0, accepted = 0, improved = 0;

    static bool between(ll a, ll x, ll b) { return a < b ? (a < x && x < b) : (x > a || x < b); }

    void respace(Rep &r) {
        // B: present non-K by bk
        std::vector<int> bl;
        for (int v = 0; v < r.n; v++) if (r.present[v] && !r.inK[v] && !r.isA2A[v]) bl.push_back(v);
        std::sort(bl.begin(), bl.end(), [&](int x, int y) { return r.bk[x] != r.bk[y] ? r.bk[x] < r.bk[y] : x < y; });
        for (size_t i = 0; i < bl.size(); i++) r.bk[bl[i]] = KEY_BASE + (ll)i * KEY_GAP;
        // A: init entries of all values + targets of present non-K
        struct E { ll k; int v; int kind; };
        std::vector<E> al;
        for (int v = 0; v < r.n; v++) {
            bool both = r.present[v] && r.inK[v];
            al.push_back({r.ik[v], v, both ? 2 : 0});
            if (r.present[v] && !r.inK[v]) al.push_back({r.tk[v], v, 1});
        }
        std::sort(al.begin(), al.end(), [](const E &x, const E &y) { return x.k != y.k ? x.k < y.k : (x.v != y.v ? x.v < y.v : x.kind < y.kind); });
        for (size_t i = 0; i < al.size(); i++) {
            ll k = KEY_BASE + (ll)i * KEY_GAP;
            if (al[i].kind == 0) r.ik[al[i].v] = k;
            else if (al[i].kind == 1) r.tk[al[i].v] = k;
            else r.ik[al[i].v] = r.tk[al[i].v] = k;
        }
    }

    void remove_value(Rep &r, int v) {
        if (!r.inK[v]) {
            r.q.erase(std::remove_if(r.q.begin(), r.q.end(), [&](const Query &x) { return x.v == v; }), r.q.end());
        }
        r.present[v] = 0;
        r.fix_minv();
    }

    // Reinsert v into r; returns delta.
    int reinsert(Rep &r, int v) {
        for (int attempt = 0; attempt < 2; attempt++) {
            // ---- B candidates ----
            std::vector<int> bl;
            for (int u = 0; u < r.n; u++) if (r.present[u] && !r.inK[u] && !r.isA2A[u]) bl.push_back(u);
            std::sort(bl.begin(), bl.end(), [&](int x, int y) { return r.bk[x] != r.bk[y] ? r.bk[x] < r.bk[y] : x < y; });
            int m = bl.size();
            std::vector<int> pos(r.n, -1);
            for (int k = 0; k < m; k++) pos[bl[k]] = k;
            std::vector<char> slotOn(std::max(m, 1), 0);
            if (m > 0) {
                if (P.ball) std::fill(slotOn.begin(), slotOn.end(), 1);
                else {
                    for (int d = 1; d <= P.bwin; d++)
                        for (int sg = -1; sg <= 1; sg += 2) {
                            int w = v + sg * d;
                            if (w < 0 || w >= r.n || pos[w] < 0) continue;
                            slotOn[pos[w]] = 1;                 // after w
                            slotOn[(pos[w] - 1 + m) % m] = 1;   // before w
                        }
                    // a couple of random slots for diversity
                    for (int k = 0; k < P.brand; k++) slotOn[rng.below(m)] = 1;
                    if (r.bk[v] != 0) {  // slots around the previous BOrder position of v
                        int p0 = (int)(std::lower_bound(bl.begin(), bl.end(), v, [&](int x, int) { return r.bk[x] < r.bk[v]; }) - bl.begin());
                        for (int d = -2; d <= 1; d++) slotOn[((p0 + d) % m + m) % m] = 1;
                    }
                }
            } else slotOn[0] = 1;
            std::vector<ll> Eb; bool bad = false;
            for (int s = 0; s < (int)slotOn.size(); s++) if (slotOn[s]) {
                ll key;
                if (m == 0) key = KEY_BASE;
                else if (s < m - 1) { ll a = r.bk[bl[s]], b = r.bk[bl[s + 1]]; if (b - a < 2) { bad = true; break; } key = a + (b - a) / 2; }
                else key = r.bk[bl[m - 1]] + KEY_GAP;
                Eb.push_back(key * 1024 + v);
            }
            // ---- A target candidates ----
            struct E { ll k; int v; };
            std::vector<E> al;
            for (int u = 0; u < r.n; u++) {
                al.push_back({r.ik[u], u});  // init entries of all values (absent ones keep their slot)
                if (r.present[u] && !r.inK[u]) al.push_back({r.tk[u], u});
            }
            std::sort(al.begin(), al.end(), [](const E &x, const E &y) { return x.k != y.k ? x.k < y.k : x.v < y.v; });
            int L = al.size();
            int pv = -1, sv = -1;
            for (int d = 1; d < r.n; d++) { int w = (v - d + r.n) % r.n; if (r.present[w] && !r.ph[w]) { pv = w; break; } }
            for (int d = 1; d < r.n; d++) { int w = (v + d) % r.n; if (r.present[w] && !r.ph[w]) { sv = w; break; } }
            std::vector<ll> Et;
            bool keepOK = false;
            if (pv < 0) {
                Et.push_back((r.ik[v] + 1) * 1024 + v);  // alone
                keepOK = true;
            } else {
                auto idx = [&](ll key) { for (int k = 0; k < L; k++) if (al[k].k == key) return k; return -1; };
                int ip = idx(r.tk[pv]), is = idx(r.tk[sv]);
                int cnt = ((is - ip) % L + L) % L;
                if (pv == sv) cnt = L;  // only one other value: every slot
                for (int c = 0; c < cnt; c++) {
                    int k = (ip + c) % L;
                    ll key;
                    if (k < L - 1) { ll a = al[k].k, b = al[k + 1].k; if (b - a < 2) { bad = true; break; } key = a + (b - a) / 2; }
                    else key = al[L - 1].k + KEY_GAP;
                    Et.push_back(key * 1024 + v);
                }
                keepOK = pv == sv ? true : between(r.tk[pv], r.ik[v], r.tk[sv]);
            }
            if (bad) { respace(r); continue; }
            ll minik = LLONG_MAX; for (int u = 0; u < r.n; u++) if (r.present[u]) minik = std::min(minik, r.ik[u]);
            bool vBeforeFirst = r.ik[v] < minik;
            if (P.twophase) {
                int sp = 0; while (sp < (int)r.q.size() && r.q[sp].t == 0) sp++;
                R.iHi = sp; R.jLo = sp;
            }
            Reinserter::Best b;
            if (P.sweep && m >= 2) {
                // all-slot sweep with provisional target, then exact refinement of the top slots
                // upper bound from value-neighbour slots (exact DP), also prepares per-target arrays
                std::vector<ll> Eub;
                for (int d = 1; d <= 2; d++) for (int sg = -1; sg <= 1; sg += 2) {
                    int w = v + sg * d;
                    if (w < 0 || w >= r.n || pos[w] < 0) continue;
                    for (int q : {pos[w], (pos[w] - 1 + m) % m}) {
                        ll key = q < m - 1 ? r.bk[bl[q]] + (r.bk[bl[q + 1]] - r.bk[bl[q]]) / 2 : r.bk[bl[m - 1]] + KEY_GAP;
                        Eub.push_back(key * 1024 + v);
                    }
                }
                std::sort(Eub.begin(), Eub.end()); Eub.erase(std::unique(Eub.begin(), Eub.end()), Eub.end());
                Reinserter::Best ub = R.solve(r, v, Eub, Et, keepOK, vBeforeFirst);
                R.ubInit = ub.delta == INT_MAX ? INT_MAX : ub.delta + 1;
                int ti0 = std::min(std::max(P.t0, 0), (int)Et.size() - 1);
                R.sweep(r, v, m, pos, ti0, vBeforeFirst);
                std::vector<std::pair<int, int>> sc;
                for (int q = 0; q < m; q++) sc.push_back({R.slotBest[q], q});
                int k = std::min(m, P.topb);
                std::partial_sort(sc.begin(), sc.begin() + k, sc.end());
                std::vector<int> sel; for (int q = 0; q < k; q++) sel.push_back(sc[q].second);
                std::sort(sel.begin(), sel.end());
                std::vector<ll> Eb2;
                for (int q : sel) {
                    ll key;
                    if (q < m - 1) { ll a = r.bk[bl[q]], bb = r.bk[bl[q + 1]]; key = a + (bb - a) / 2; }
                    else key = r.bk[bl[m - 1]] + KEY_GAP;
                    Eb2.push_back(key * 1024 + v);
                }
                R.ubInit = INT_MAX;
                int tm = R.tmode; R.tmode = 0;
                b = R.solve(r, v, Eb2, Et, keepOK, vBeforeFirst);
                R.tmode = tm;
                Eb = Eb2;
                if (ub.delta < b.delta) { b = ub; Eb = Eub; }
            } else b = R.solve(r, v, Eb, Et, keepOK, vBeforeFirst);
            lastInfo = b; lastQ = r.q.size(); lastVBF = vBeforeFirst;
            r.present[v] = 1;
            r.isA2A[v] = 0;
            if (b.keep) { r.inK[v] = 1; r.tk[v] = r.ik[v]; }
            else if (b.a2a) {
                r.inK[v] = 0; r.isA2A[v] = 1; r.tk[v] = Et[b.ti] >> 10;
                r.q.insert(r.q.begin() + b.i, {(int16_t)v, 2});
            } else {
                r.inK[v] = 0;
                r.bk[v] = Eb[b.bi] >> 10; r.tk[v] = Et[b.ti] >> 10;
                r.q.insert(r.q.begin() + b.i, {(int16_t)v, 0});
                r.q.insert(r.q.begin() + b.j + 1, {(int16_t)v, 1});
            }
            r.fix_minv();
            return b.delta;
        }
        fprintf(stderr, "reinsert failed\n"); exit(4);
    }

    void compute_vcost(const Rep &r) {
        vcost.assign(r.n, 0);
        St s = initial_state(r);
        for (auto qu : r.q) vcost[qu.v] += do_step(s, r, qu);
    }

    std::vector<int> pick_destroy(const Rep &r) {
        std::vector<char> mark(r.n, 0);
        std::vector<int> D;
        int mode = rng.below(10);
        int Q = r.q.size();
        if (mode < 6 && Q > 0) {
            int L = 1 + rng.below(P.seglen);
            int st = rng.below(Q);
            for (int k = st; k < std::min(Q, st + L); k++) { int v = r.q[k].v; if (!mark[v]) { mark[v] = 1; D.push_back(v); } }
        } else if (mode < 8) {
            int w = 1 + rng.below(P.vwin), x = rng.below(r.n);
            for (int k = 0; k < w; k++) { int v = (x + k) % r.n; if (!mark[v]) { mark[v] = 1; D.push_back(v); } }
        } else {
            int c = 1 + rng.below(P.rmax);
            for (int k = 0; k < c; k++) { int v = rng.below(r.n); if (!mark[v]) { mark[v] = 1; D.push_back(v); } }
        }
        if ((int)D.size() >= r.n - 1) D.resize(std::max(1, r.n - 2));
        return D;
    }

    void run(const Rep &init, double tlimit) {
        cur = init; curCost = full_cost(cur); best = cur; bestCost = curCost;
        std::vector<int> hist(std::max(1, P.lahc), curCost);
        double t0 = now_sec();
        for (iters = 0;; iters++) {
            if (P.iters >= 0 ? iters >= P.iters : ((iters & 7) == 0 && now_sec() - t0 > tlimit)) break;
            Rep w = cur;
            compute_vcost(w);
            std::vector<int> D = pick_destroy(w);
            std::vector<int> dc(D.size());
            for (size_t k = 0; k < D.size(); k++) dc[k] = vcost[D[k]];
            std::vector<int> ord(D.size());
            for (size_t k = 0; k < D.size(); k++) ord[k] = k;
            std::sort(ord.begin(), ord.end(), [&](int x, int y) { return dc[x] != dc[y] ? dc[x] > dc[y] : x < y; });
            bool qd = rng.uni() < P.pqd;
            for (int v : D) {
                if (qd) {
                    w.q.erase(std::remove_if(w.q.begin(), w.q.end(), [&](const Query &x) { return x.v == v; }), w.q.end());
                    w.inK[v] = 1; w.isA2A[v] = 0; w.tk[v] = w.ik[v]; w.ph[v] = 1;
                } else remove_value(w, v);
            }
            int c = full_cost(w);
            for (int k : ord) {
                if (qd) { w.ph[D[k]] = 0; remove_value(w, D[k]); c = full_cost(w); }
                Rep before = w;
                int d = reinsert(w, D[k]);
                c += d;
                if (P.a2a && !P.check) c = full_cost(w);
                if (P.check) {
                    int fc = full_cost(w);
                    if (fc != c && lastInfo.keep) {
                        int v = D[k]; ll EkI = w.ik[v] * 1024 + v;
                        St s = initial_state(w);
                        for (int i = 0; i < (int)w.q.size(); i++) {
                            int real = do_step(s, w, w.q[i]);
                            int topo = (i == 0 && lastVBF) ? R.rkI[0] : -1;
                            int pred = R.aug(i, EkI, R.rkI[i], -1, 0, topo, -1);
                            if (real != pred) { auto &f = R.qi[i]; fprintf(stderr, "q%d t=%d u=%d real=%d pred=%d old=%d sa=%d sb=%d lbA=%d srcIdx=%d dstRaw=%d rkI=%d ancA=%lld EkI=%lld src=%lld dst=%lld\n", i, f.t, w.q[i].v, real, pred, f.c, f.sa, f.sb, f.lbA, f.srcIdx, f.dstRaw, R.rkI[i], f.ancA>>10, EkI>>10, f.src>>10, f.dst>>10); }
                        }
                    }
                    if (fc != c && !lastInfo.keep && !lastInfo.a2a && !P.a2a) { R.diag(w, D[k], lastInfo, lastVBF); fprintf(stderr, "c_before=%d full_before=%d delta=%d\n", c - d, full_cost(before), d); }
                    if (fc != c && P.a2a) { nMis++; c = fc; }
                    if (fc != c) { fprintf(stderr, "DELTA MISMATCH v=%d pred=%d real=%d keep=%d i=%d j=%d Q=%d vbf=%d minv=%d\n", D[k], c, fc, lastInfo.keep, lastInfo.i, lastInfo.j, lastQ, lastVBF, w.minv); exit(5); }
                }
            }
            int h = hist[iters % hist.size()];
            if (c <= curCost || c <= h) {
                if (c < curCost) improved++;
                cur = std::move(w); curCost = c; accepted++;
                if (c < bestCost) { bestCost = c; best = cur; }
            }
            hist[iters % hist.size()] = curCost;
        }
    }
};

int main(int argc, char **argv) {
    Params P;
    // size-dependent defaults (overridable): tuned on train/dev sets, see EXPERIMENTS.md
    P.sweep = 1; P.pqd = 0.5; P.time = -1;
    int first = 1;
    bool setLarge = false, setPeep = false;
    for (; first < argc && argv[first][0] == '-' && argv[first][1] == '-'; first++) {
        std::string a = argv[first];
        auto val = [&](const char *k) -> const char * { size_t L = strlen(k); return a.compare(0, L, k) == 0 ? a.c_str() + L : nullptr; };
        if (auto x = val("--time=")) P.time = atof(x);
        else if (auto x = val("--iters=")) P.iters = atoll(x);
        else if (auto x = val("--seed=")) P.seed = strtoull(x, 0, 10);
        else if (auto x = val("--chunks=")) P.chunks = atoi(x);
        else if (auto x = val("--bwin=")) P.bwin = atoi(x);
        else if (auto x = val("--ball=")) P.ball = atoi(x);
        else if (auto x = val("--tmode=")) P.tmode = atoi(x);
        else if (auto x = val("--twophase=")) P.twophase = atoi(x);
        else if (auto x = val("--brand=")) P.brand = atoi(x);
        else if (auto x = val("--sweep=")) P.sweep = atoi(x);
        else if (auto x = val("--a2a=")) P.a2a = atoi(x);
        else if (auto x = val("--pqd=")) P.pqd = atof(x);
        else if (auto x = val("--restarts=")) P.restarts = atoi(x);
        else if (auto x = val("--rfrac=")) P.rfrac = atof(x);
        else if (auto x = val("--kinit=")) P.kinit = atoi(x);
        else if (auto x = val("--peep=")) { P.peepK.clear(); for (const char *q = x; *q;) { P.peepK.push_back(atoi(q)); while (*q && *q != ',') q++; if (*q) q++; } }
        else if (auto x = val("--a2amax=")) P.a2amax = atoi(x);
        else if (auto x = val("--ginit=")) P.ginit = atoi(x);
        else if (auto x = val("--bw=")) P.bw = atoi(x);
        else if (auto x = val("--bc=")) P.bc = atoi(x);
        else if (auto x = val("--nlis=")) P.nlis = atoi(x);
        else if (auto x = val("--bws=")) P.bws = atoi(x);
        else if (auto x = val("--ktop=")) P.ktop = atoi(x);
        else if (auto x = val("--t0=")) P.t0 = atoi(x);
        else if (auto x = val("--topb=")) P.topb = atoi(x);
        else if (auto x = val("--seglen=")) P.seglen = atoi(x);
        else if (auto x = val("--rmax=")) P.rmax = atoi(x);
        else if (auto x = val("--vwin=")) P.vwin = atoi(x);
        else if (auto x = val("--lahc=")) P.lahc = atoi(x);
        else if (auto x = val("--verbose=")) P.verbose = atoi(x);
        else if (auto x = val("--check=")) P.check = atoi(x);
        else { fprintf(stderr, "Error\n"); return 1; }
        if (a.rfind("--peep", 0) == 0) setPeep = true;
        if (a.rfind("--seglen", 0) == 0 || a.rfind("--nlis", 0) == 0 || a.rfind("--bw", 0) == 0) setLarge = true;
    }
    std::vector<long long> vals;
    if (!parse_args(argc, argv, first, vals)) { fprintf(stderr, "Error\n"); return 1; }
    int n = vals.size();
    if (n == 0) return 0;
    std::vector<int> a0 = to_ranks(vals);
    bool sorted = true; for (int i = 0; i < n; i++) if (a0[i] != i) sorted = false;
    if (sorted) return 0;
    if (n <= 8) { print_ops(exact::solve(a0)); return 0; }
    if (n > 1000) {  // representation encodes values in 10 bits: constructive fallback
        Rng rr(1); auto K = circular_lis(a0, rr);
        print_ops(chunk_construct(a0, K, std::max(2, n / 80), 0)); return 0;
    }
    if (P.time < 0) P.time = n <= 150 ? 10.0 : (n <= 300 ? 20.0 : 40.0);
    if (n > 300 && !setLarge) { P.seglen = 16; P.rmax = 9; P.vwin = 9; P.nlis = 30; P.ktop = 3; P.bw = 100; P.bws = 0; }
    if (!setPeep) P.peepK = n <= 150 ? std::vector<int>{4, 6, 8} : std::vector<int>{4, 6};
    double t0g = now_sec();
    Rng rng(P.seed * 0x9E3779B97F4A7C15ull + 12345);
    int m = P.chunks > 0 ? P.chunks : (n <= 150 ? 2 : (n <= 300 ? 4 : 6));
    int lis_len;
    std::vector<Op> base;
    std::vector<std::vector<Op>> inits;
    {
        // candidate K sets: LIS from many start positions (distinct), screened with a small beam
        std::vector<std::vector<char>> Ks;
        int starts = std::min(n, P.nlis);
        for (int t = 0; t < starts; t++) {
            int st = (int)((long long)t * n / starts);
            auto K = lis_from(a0, st);
            if (std::find(Ks.begin(), Ks.end(), K) == Ks.end()) Ks.push_back(K);
        }
        std::vector<std::pair<size_t, int>> sc;
        for (int k = 0; k < (int)Ks.size(); k++) sc.push_back({chunk_construct(a0, Ks[k], m, P.bws, P.bc).size(), k});
        std::sort(sc.begin(), sc.end());
        for (int k = 0; k < std::min((int)sc.size(), std::max(P.ktop, P.kinit)); k++) {
            std::vector<Op> o = chunk_construct(a0, Ks[sc[k].second], m, P.bw, P.bc);
            inits.push_back(o);
            if (base.empty() || o.size() < base.size()) base = o;
        }
        lis_len = 0; for (char c : Ks[sc[0].second]) lis_len += c;
    }
    Rep rep;
    if (!verify_ops(a0, base) || !convert(a0, base, rep)) { fprintf(stderr, "internal: constructor\n"); return 2; }
    if (P.verbose) fprintf(stderr, "init=%zu valid=%d t=%.2f\n", base.size(), (int)verify_ops(a0, base), now_sec() - t0g);
    IG ig; ig.P = P; ig.rng = Rng(P.seed + 99);
    ig.R.tmode = P.tmode; ig.R.t0 = P.t0; ig.R.topb = P.topb; ig.R.allowA2A = P.a2a; ig.R.a2aMax = P.a2amax;
    if (P.ginit) {
        Rep g = rep;
        std::vector<int> order;
        for (int v = 0; v < n; v++) if (!g.inK[v]) order.push_back(v);
        for (int v : order) ig.remove_value(g, v);
        if (P.ginit == 1) std::sort(order.begin(), order.end(), [](int a, int b) { return a > b; });
        else for (int i = (int)order.size() - 1; i > 0; i--) std::swap(order[i], order[ig.rng.below(i + 1)]);
        for (int v : order) { g.bk[v] = 0; ig.reinsert(g, v); }
        int gc = full_cost(g);
        if (P.verbose) fprintf(stderr, "ginit=%d t=%.2f\n", gc, now_sec() - t0g);
        if (gc < full_cost(rep)) rep = g;
    }
    double t0 = now_sec();
    // budget: total --time covers init + IG + peephole (reserve for peephole)
    P.time = std::max(0.5, P.time - (t0 - t0g) - (n <= 150 ? 2.0 : 4.0));
    ig.P.time = P.time;
    if (P.restarts > 1 && P.iters < 0) {
        // racing: R independent runs for frac*T, then continue the best for the rest
        double slice = P.time * P.rfrac / P.restarts;
        Rep bestRep = rep; int bestC = INT_MAX;
        for (int k = 0; k < P.restarts; k++) {
            IG g; g.P = P; g.rng = Rng(P.seed * 1000 + 99 + k * 7919); g.R = ig.R;
            Rep start = rep;
            if (P.kinit > 1 && !inits.empty()) {
                Rep alt; const std::vector<Op> &o = inits[k % std::min((int)inits.size(), P.kinit)];
                if (convert(a0, o, alt)) start = alt;
            }
            g.run(start, slice);
            if (g.bestCost < bestC) { bestC = g.bestCost; bestRep = g.best; }
            ig.iters += g.iters;
        }
        long long it0 = ig.iters;
        ig.run(bestRep, std::max(0.0, P.time - (now_sec() - t0)));
        ig.iters += it0;
    } else
    ig.run(rep, P.time);
    std::vector<Op> ops = emit(ig.best);
    if ((int)ops.size() != ig.bestCost || !verify_ops(a0, ops)) {
        fprintf(stderr, "internal: emit mismatch %zu vs %d\n", ops.size(), ig.bestCost);
        ops = base;
    }
    int psaved = 0;
    for (int K : P.peepK) { if (K >= 2) psaved += peep::optimize(a0, ops, K); }
    if (!verify_ops(a0, ops)) { fprintf(stderr, "internal: peephole broke solution\n"); ops = emit(ig.best); }
    if (P.verbose) fprintf(stderr, "peephole saved %d\n", psaved);
    if (P.verbose)
        fprintf(stderr, "n=%d lis=%d constructor=%zu final=%zu time=%.2f iters=%lld acc=%lld impr=%lld mis=%lld a2a=%d\n", n, lis_len,
                base.size(), ops.size(), now_sec() - t0, ig.iters, ig.accepted, ig.improved, ig.nMis,
                (int)std::count(ig.best.isA2A.begin(), ig.best.isA2A.end(), 1));
    print_ops(ops);
    return 0;
}
