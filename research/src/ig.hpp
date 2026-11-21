// Iterated Greedy over the query representation with exact best re-insertion of one value.
//
// Re-insertion of value v (currently absent) chooses
//   i : A2B(v) inserted before old query i
//   j : B2A(v) inserted before old query j (j >= i; if j == i, A2B then B2A both before query i)
//   b : BOrder slot (candidate set), t : AOrder target slot (candidate set)
// or "keep v in A" (K option). For fixed (b,t) the best (i,j) is found by an O(Q) DP using
// the fact that v's presence changes every other query's rotation by a computable amount
// (each query's cost is re-evaluated exactly with v inserted, through anchors for the tops).
#pragma once
#include "rep.hpp"

struct QInfo {
    int sa, sb;          // stack sizes before the query
    ll ancA, ancB;       // top anchors (encoded): top = first element >= anchor (cyclic)
    int lbA, lbB;        // lower_bound of anchors in A,B
    ll src, dst;         // encoded source key (element moved) and destination key
    int srcIdx, dstRaw;  // index of source element in its stack, raw lower_bound of dst in dest
    int t;               // 0 A2B 1 B2A
    int c;               // cost
};

struct Reinserter {
    // scratch buffers
    std::vector<QInfo> qi;
    std::vector<int> rkI;                      // rank of kI in A_i
    std::vector<std::vector<int>> rb, rt;      // ranks of candidate keys
    std::vector<int> P, G, SUF, Aarr;
    int Q = 0, finSa = 0, finTa = 0, finOld = 0;
    ll finAnc = 0; int finLb = 0;

    static inline int md(int x, int m) { return m > 0 ? ((x % m) + m) % m : 0; }
    // fast cyclic reduction for x in (-m, 2m)
    static inline int mq(int x, int m) { x = x < 0 ? x + m : x; x = x >= m ? x - m : x; return m > 0 ? x : 0; }

    // Pass over the current solution recording per-query information and candidate ranks.
    void pass(const Rep &r, ll EkI, const std::vector<ll> &Eb, const std::vector<ll> &Et) {
        Q = r.q.size();
        bool ebSorted = std::is_sorted(Eb.begin(), Eb.end());
        qi.resize(Q);
        rkI.resize(Q + 1);
        rb.assign(Eb.size(), std::vector<int>(Q + 1));
        rt.assign(Et.size(), std::vector<int>(Q + 1));
        St s = initial_state(r);
        // initial anchors: exact top element
        ll ancA = s.A.empty() ? 0 : s.A[s.ta], ancB = 0;
        for (int i = 0; i <= Q; i++) {
            int sa = s.A.size(), sb = s.B.size();
            rkI[i] = (int)(std::lower_bound(s.A.begin(), s.A.end(), EkI) - s.A.begin());
            if (ebSorted) {
                size_t p = 0, nbk = Eb.size(), L = s.B.size();
                for (size_t k = 0; k < nbk; k++) { while (p < L && s.B[p] < Eb[k]) p++; rb[k][i] = (int)p; }
            } else
                for (size_t k = 0; k < Eb.size(); k++) rb[k][i] = (int)(std::lower_bound(s.B.begin(), s.B.end(), Eb[k]) - s.B.begin());
            for (size_t k = 0; k < Et.size(); k++) rt[k][i] = (int)(std::lower_bound(s.A.begin(), s.A.end(), Et[k]) - s.A.begin());
            if (i == Q) {
                finSa = sa; finTa = s.ta; finAnc = ancA;
                finLb = (int)(std::lower_bound(s.A.begin(), s.A.end(), ancA) - s.A.begin());
                finOld = final_rot(s, r);
                minIdx = r.minv < r.n ? (int)(std::lower_bound(s.A.begin(), s.A.end(), r.tk[r.minv] * 1024 + r.minv) - s.A.begin()) : 0;
                break;
            }
            QInfo &f = qi[i];
            f.sa = sa; f.sb = sb; f.ancA = ancA; f.ancB = ancB;
            f.lbA = (int)(std::lower_bound(s.A.begin(), s.A.end(), ancA) - s.A.begin());
            f.lbB = (int)(std::lower_bound(s.B.begin(), s.B.end(), ancB) - s.B.begin());
            Query qu = r.q[i]; int u = qu.v; f.t = qu.t;
            if (qu.t == 0) {
                f.src = r.ik[u] * 1024 + u; f.dst = r.bk[u] * 1024 + u;
                f.srcIdx = (int)(std::lower_bound(s.A.begin(), s.A.end(), f.src) - s.A.begin());
                f.dstRaw = (int)(std::lower_bound(s.B.begin(), s.B.end(), f.dst) - s.B.begin());
            } else {
                f.src = r.bk[u] * 1024 + u; f.dst = r.tk[u] * 1024 + u;
                f.srcIdx = (int)(std::lower_bound(s.B.begin(), s.B.end(), f.src) - s.B.begin());
                f.dstRaw = (int)(std::lower_bound(s.A.begin(), s.A.end(), f.dst) - s.A.begin());
            }
            f.c = do_step(s, r, qu);
            // anchors for the next query
            if (qu.t == 0) { ancA = f.src; ancB = f.dst; }   // A: successor of popped; B: exact pushed
            else { ancA = f.dst; ancB = f.src; }
        }
    }

    // cost of old query i with an optional extra element in A (xa, rank ra) and/or B (xb, rank rbk).
    // topA/topB overrides: -1 = use anchor.
    inline int aug(int i, ll xa, int ra, ll xb, int rbk, int topAo = -1, int topBo = -1) const {
        const QInfo &f = qi[i];
        int sa = f.sa + (xa >= 0), sb = f.sb + (xb >= 0);
        int ta = topAo >= 0 ? topAo : md(f.lbA + (xa >= 0 && xa < f.ancA), sa);
        int tb = topBo >= 0 ? topBo : md(f.lbB + (xb >= 0 && xb < f.ancB), sb);
        (void)ra; (void)rbk;
        int fa, fb;
        if (f.t == 0) {
            int si = f.srcIdx + (xa >= 0 && xa < f.src);
            fa = md(si - ta, sa);
            int sl = sb ? md(f.dstRaw + (xb >= 0 && xb < f.dst), sb) : 0;
            fb = sb ? md(sl - tb, sb) : 0;
        } else {
            int si = f.srcIdx + (xb >= 0 && xb < f.src);
            fb = md(si - tb, sb);
            int sl = sa ? md(f.dstRaw + (xa >= 0 && xa < f.dst), sa) : 0;
            fa = sa ? md(sl - ta, sa) : 0;
        }
        return 1 + comb_cost(fa, sa, fb, sb);
    }

    struct Best { int delta; int i, j, bi, ti; bool keep; };

    // r: current solution (v absent). Candidate encoded keys for B (Eb) and A targets (Et).
    int tmode = 1, t0 = 0, topb = 3;
    bool tOnly = false;
    int iHi = INT_MAX, jLo = 0;  // constraints: A2B before old query i <= iHi, B2A before old query j >= jLo
    Best solve(const Rep &r, int v, const std::vector<ll> &Eb, const std::vector<ll> &Et, bool keepOK, bool vBeforeFirst) {
        ll EkI = r.ik[v] * 1024 + v;
        pass(r, EkI, Eb, Et);
        Best best{INT_MAX, -1, -1, -1, -1, false};
        // prefix effect of v in A at kI
        P.assign(Q + 1, 0);
        for (int i = 0; i < Q; i++) {
            int topo = -1;
            if (i == 0 && vBeforeFirst) topo = rkI[0];
            int c2 = aug(i, EkI, rkI[i], -1, 0, topo, -1);
            P[i + 1] = P[i] + c2 - qi[i].c;
        }
        // final rotation with an extra element x (rank rx) in A, top override
        auto fin = [&](ll x, int rx, int topo) {
            int sa = finSa + 1;
            int ta = topo >= 0 ? topo : md(finLb + (x < finAnc), sa);
            int tgt;
            if (v < r.minv) tgt = rx;
            else { ll em = r.tk[r.minv] * 1024 + r.minv; tgt = minIdx + (x < em); }
            int f = md(tgt - ta, sa);
            return std::min(f, sa - f) - finOld;
        };
        if (keepOK) {
            int topo = (Q == 0 && vBeforeFirst) ? rkI[0] : -1;
            int d = P[Q] + fin(EkI, rkI[Q], topo);
            if (d < best.delta) best = {d, -1, -1, -1, -1, true};
        }
        int nb = Eb.size(), nt = Et.size();
        // ---- per-t precomputation ----
        tSUF.resize(nt); tFaC2.resize(nt); tFaS2.resize(nt); tFaCC.resize(nt); tFinA.resize(nt); tFinB.resize(nt);
        for (int ti = 0; ti < nt; ti++) {
            ll et = Et[ti];
            const std::vector<int> &T = rt[ti];
            std::vector<int> &SU = tSUF[ti]; SU.assign(Q + 2, 0);
            for (int q = Q - 1; q >= 0; q--) SU[q] = SU[q + 1] + aug(q, et, T[q], -1, 0) - qi[q].c;
            tFinA[ti] = fin(et, T[Q], -1); tFinB[ti] = fin(et, T[Q], T[Q]);
            std::vector<int> &F2 = tFaC2[ti], &FS = tFaS2[ti], &FC = tFaCC[ti];
            F2.resize(Q + 1); FS.resize(Q + 1); FC.resize(Q + 1);
            for (int j = 0; j <= Q; j++) {
                int sa = j < Q ? qi[j].sa : finSa;
                int ta = j < Q ? md(qi[j].lbA, sa) : md(finLb, finSa);
                F2[j] = sa ? md(md(T[j], sa) - ta, sa) : 0;
                FS[j] = j < Q ? augA(j, et, T[j]) : 0;
                int tA = sa ? md(rkI[j], sa) : 0;
                FC[j] = sa ? md(md(T[j], sa) - tA, sa) : 0;
            }
        }
        if (tOnly) return best;
        G.resize(Q + 1); Aarr.resize(Q + 1);
        bFbC2.resize(Q + 1); bFbS2.resize(Q + 1); bC1ii.resize(Q + 1);
        auto prep_b = [&](int bi) {
            ll eb = Eb[bi];
            const std::vector<int> &R = rb[bi];
            G[0] = 0;
            for (int q = 0; q < Q; q++) G[q + 1] = G[q] + aug(q, -1, 0, eb, R[q]) - qi[q].c;
            for (int i = 0; i <= Q; i++) {
                int sa = i < Q ? qi[i].sa : finSa, sb = i < Q ? qi[i].sb : 0;
                int sa1 = sa + 1;
                int ta1;
                if (i == 0 && vBeforeFirst) ta1 = rkI[0];
                else ta1 = i < Q ? md(qi[i].lbA + (EkI < qi[i].ancA), sa1) : md(finLb + (EkI < finAnc), sa1);
                int fa = md(rkI[i] - ta1, sa1);
                int fb = (i < Q && sb) ? md(md(R[i], sb) - md(qi[i].lbB, sb), sb) : 0;
                int c1 = 1 + comb_cost(fa, sa1, fb, sb);
                bC1ii[i] = c1;
                if (i < Q) {
                    const QInfo &f = qi[i];
                    int s1 = aug(i, -1, 0, eb, R[i], f.sa ? md(rkI[i], f.sa) : 0, R[i]);
                    Aarr[i] = P[i] + c1 + s1 - f.c - G[i + 1];
                    int sb1 = f.sb + 1;
                    int tb1 = md(f.lbB + (eb < f.ancB), sb1);
                    bFbC2[i] = md(R[i] - tb1, sb1);
                    bFbS2[i] = augB(i, -1, f.sb ? md(R[i], f.sb) : 0);
                } else { bFbC2[i] = 0; bFbS2[i] = 0; }
            }
        };
        auto dp_bt = [&](int bi, int ti, Best &bst) {
            const std::vector<int> &SU = tSUF[ti], &F2 = tFaC2[ti], &FS = tFaS2[ti], &FC = tFaCC[ti];
            int finA = tFinA[ti], finB = tFinB[ti];
            int pm = INT_MAX, pmi = -1;
            for (int j = 0; j <= Q; j++) {
                if (j >= 1 && j - 1 <= iHi && Aarr[j - 1] < pm) { pm = Aarr[j - 1]; pmi = j - 1; }
                if (j < jLo) continue;
                int tail, c2, s2 = 0;
                if (j < Q) {
                    const QInfo &f = qi[j];
                    c2 = 1 + comb_cost(F2[j], f.sa, bFbC2[j], f.sb + 1);
                    s2 = 1 + comb_cost(FS[j], f.sa + 1, bFbS2[j], f.sb);
                    tail = s2 - f.c + SU[j + 1] + finA;
                } else {
                    c2 = 1 + comb_cost(F2[Q], finSa, 0, 1);
                    tail = finB;
                }
                int Bj = G[j] + c2 + tail;
                if (pmi >= 0 && pm + Bj < bst.delta) bst = {pm + Bj, pmi, j, bi, ti, false};
                int sa = j < Q ? qi[j].sa : finSa, sb = j < Q ? qi[j].sb : 0;
                int cc = 1 + comb_cost(FC[j], sa, 0, sb + 1);
                int d = P[j] + bC1ii[j] + cc + tail;
                if (j <= iHi && d < bst.delta) bst = {d, j, j, bi, ti, false};
            }
        };
        if (tmode == 0 || nt == 1) {
            for (int bi = 0; bi < nb; bi++) { prep_b(bi); for (int ti = 0; ti < nt; ti++) dp_bt(bi, ti, best); }
        } else {
            int tp = std::min(std::max(t0, 0), nt - 1);
            std::vector<std::pair<int, int>> cand;  // (delta, bi)
            for (int bi = 0; bi < nb; bi++) {
                prep_b(bi);
                Best bb{INT_MAX, -1, -1, -1, -1, false};
                dp_bt(bi, tp, bb);
                cand.push_back({bb.delta, bi});
                if (bb.delta < best.delta) best = bb;
            }
            int k = std::min((int)cand.size(), topb);
            std::partial_sort(cand.begin(), cand.begin() + k, cand.end());
            for (int c = 0; c < k; c++) {
                prep_b(cand[c].second);
                for (int ti = 0; ti < nt; ti++) dp_bt(cand[c].second, ti, best);
            }
        }
        return best;
    }
    // ---- all-slot sweep: B slots s = 0..m-1 (slot s lies between bl[s] and bl[s+1]); posOf[value] = bl index.
    // Evaluates every slot with provisional target ti0, returns per-slot best deltas in slotBest.
    std::vector<int> swR, swG, swPm, swPmi, slotBest, slotI, slotJ;
    Best sweep(const Rep &r, int v, int m, const std::vector<int> &posOf, int ti0, bool vBeforeFirst) {
        ll EkI = r.ik[v] * 1024 + v;
        const std::vector<int> &SU = tSUF[ti0], &F2 = tFaC2[ti0], &FS = tFaS2[ti0], &FC = tFaCC[ti0];
        int finA = tFinA[ti0], finB = tFinB[ti0];
        swR.assign(m, 0); swG.assign(m, 0); swPm.assign(m, INT_MAX); swPmi.assign(m, -1);
        slotBest.assign(m, INT_MAX); slotI.assign(m, -1); slotJ.assign(m, -1);
        Best best{INT_MAX, -1, -1, -1, -1, false};
        auto pv = [&](ll key) { return posOf[key & 1023]; };
        for (int x = 0; x <= Q; x++) {
            int sa = x < Q ? qi[x].sa : finSa, sb = x < Q ? qi[x].sb : 0;
            int sa1 = sa + 1;
            int ta1 = (x == 0 && vBeforeFirst) ? rkI[0] : (x < Q ? md(qi[x].lbA + (EkI < qi[x].ancA), sa1) : md(finLb + (EkI < finAnc), sa1));
            int fa1 = md(rkI[x] - ta1, sa1);
            int tbOld = x < Q ? md(qi[x].lbB, sb) : 0;
            int tA = sa ? md(rkI[x], sa) : 0;
            int cc = 1 + comb_cost(FC[x], sa, 0, sb + 1);
            if (x == Q) {
                int c2 = 1 + comb_cost(F2[Q], finSa, 0, 1);
                for (int s = 0; s < m; s++) {
                    int c1 = 1 + comb_cost(fa1, sa1, 0, 0);
                    int d = P[Q] + c1 + cc + finB;
                    if (d < slotBest[s]) { slotBest[s] = d; slotI[s] = Q; slotJ[s] = Q; }
                    if (swPmi[s] >= 0) {
                        int dd = swPm[s] + swG[s] + c2 + finB;
                        if (dd < slotBest[s]) { slotBest[s] = dd; slotI[s] = swPmi[s]; slotJ[s] = Q; }
                    }
                }
                break;
            }
            const QInfo &f = qi[x];
            int sb1 = sb + 1;
            int pA = sb ? pv(f.ancB) : -1;  // anchor position (bl index)
            int u = (f.t == 0 ? f.dst : f.src) & 1023;
            int pU = posOf[u];
            int faOld = f.t == 0 ? md(f.srcIdx - md(f.lbA, sa), sa) : (sa ? md(md(f.dstRaw, sa) - md(f.lbA, sa), sa) : 0);
            // A-part of s1: old query x with A top = element after v
            int faS1 = f.t == 0 ? md(f.srcIdx - tA, sa) : (sa ? md(md(f.dstRaw, sa) - tA, sa) : 0);
            // middle-effect e_x(s) for the 4 combinations of (s < pA, s < pU)
            int ev[2][2];
            for (int a = 0; a < 2; a++) for (int b2 = 0; b2 < 2; b2++) {
                int tb_ = sb ? md(f.lbB + a, sb1) : 0;
                int fb;
                if (f.t == 0) { int sl = md(f.dstRaw + b2, sb1); fb = md(sl - tb_, sb1); }
                else { int si = f.srcIdx + b2; fb = md(si - tb_, sb1); }
                ev[a][b2] = 1 + comb_cost(faOld, sa, fb, sb1) - f.c;
            }
            int base = SU[x + 1] + finA - f.c;
            for (int s = 0; s < m; s++) {
                int R = swR[s];
                int ltA = s < pA, ltU = s < pU;
                // ---- j = x ----
                int tb1 = sb ? mq(f.lbB + ltA, sb1) : 0;
                int fbC2 = mq(R - tb1, sb1);
                int c2 = 1 + comb_cost(F2[x], sa, fbC2, sb1);
                int fbS2 = 0;
                if (sb) { int Rm = R < sb ? R : R - sb; fbS2 = f.t == 0 ? mq(mq(f.dstRaw, sb) - Rm, sb) : mq(f.srcIdx - Rm, sb); }
                int s2 = 1 + comb_cost(FS[x], sa + 1, fbS2, sb);
                int tail = s2 + base;
                int fb1 = sb ? mq((R < sb ? R : R - sb) - tbOld, sb) : 0;
                int c1 = 1 + comb_cost(fa1, sa1, fb1, sb);
                int d = P[x] + c1 + cc + tail;
                if (d < slotBest[s]) { slotBest[s] = d; slotI[s] = x; slotJ[s] = x; }
                if (swPmi[s] >= 0) {
                    int dd = swPm[s] + swG[s] + c2 + tail;
                    if (dd < slotBest[s]) { slotBest[s] = dd; slotI[s] = swPmi[s]; slotJ[s] = x; }
                }
                // ---- A_s(x) ----
                int fbs1;
                if (f.t == 0) { int sl = mq(f.dstRaw + ltU, sb1); fbs1 = mq(sl - R, sb1); }
                else { int si = f.srcIdx + ltU; fbs1 = mq(si - R, sb1); }
                int s1 = 1 + comb_cost(faS1, sa, fbs1, sb1);
                int e = ev[ltA][ltU];
                int Ax = P[x] + c1 + s1 - f.c - (swG[s] + e);
                if (Ax < swPm[s]) { swPm[s] = Ax; swPmi[s] = x; }
                swG[s] += e;
            }
            // rank update
            if (f.t == 0) { for (int s = pU; s < m; s++) swR[s]++; }
            else { for (int s = pU; s < m; s++) swR[s]--; }
        }
        for (int s = 0; s < m; s++) if (slotBest[s] < best.delta) best = {slotBest[s], slotI[s], slotJ[s], s, ti0, false};
        return best;
    }

    int minIdx = 0;
    std::vector<std::vector<int>> tSUF, tFaC2, tFaS2, tFaCC;
    std::vector<int> tFinA, tFinB, bFbC2, bFbS2, bC1ii;

    // A-side forward distance of old query i with extra element xa in A and A top override.
    inline int augA(int i, ll xa, int topAo) const {
        const QInfo &f = qi[i];
        int sa = f.sa + 1;
        int ta = topAo;
        if (f.t == 0) { int si = f.srcIdx + (xa < f.src); return md(si - ta, sa); }
        int sl = md(f.dstRaw + (xa < f.dst), sa);
        return md(sl - ta, sa);
    }
    // B-side forward distance of old query i (no extra element in B) with B top override.
    inline int augB(int i, ll xb, int topBo) const {
        (void)xb;
        const QInfo &f = qi[i];
        int sb = f.sb;
        if (!sb) return 0;
        if (f.t == 0) { int sl = md(f.dstRaw, sb); return md(sl - topBo, sb); }
        return md(f.srcIdx - topBo, sb);
    }

    // Debug: compare per-query predicted costs against a real decode of the solution after insertion.
    void diag(const Rep &w, int v, const Best &b, bool vbf) {
        ll EkI = w.ik[v] * 1024 + v, eb = w.bk[v] * 1024 + v, et = w.tk[v] * 1024 + v;
        const std::vector<int> &R = rb[b.bi], &T = rt[b.ti];
        St s = initial_state(w);
        int nq = w.q.size(); long long sumReal = 0, oldT = finOld;
        for (int q = 0; q < Q; q++) oldT += qi[q].c;
        for (int k = 0; k < nq; k++) {
            int real = do_step(s, w, w.q[k]); sumReal += real;
            int pred = -1; const char *what = "";
            int i = b.i, j = b.j;
            if (k < i) { pred = aug(k, EkI, rkI[k], -1, 0, (k == 0 && vbf) ? rkI[0] : -1, -1); what = "pre"; }
            else if (k == i) { what = "c1";
                const QInfo &f = qi[i]; int sa1 = f.sa + 1;
                int ta1 = (i == 0 && vbf) ? rkI[0] : md(f.lbA + (EkI < f.ancA), sa1);
                int fa = md(rkI[i] - ta1, sa1);
                int fb = f.sb ? md(md(R[i], f.sb) - md(f.lbB, f.sb), f.sb) : 0;
                pred = 1 + comb_cost(fa, sa1, fb, f.sb); }
            else if (k == i + 1 && j > i) { what = "s1"; const QInfo &f = qi[i]; pred = aug(i, -1, 0, eb, R[i], f.sa ? md(rkI[i], f.sa) : 0, R[i]); }
            else if (k <= j && j > i) { what = "mid"; pred = aug(k - 1, -1, 0, eb, R[k - 1]); }
            else if (k == j + 1) { what = "c2";
                if (j < Q) { const QInfo &f = qi[j]; int sb1 = f.sb + 1; int tb1 = md(f.lbB + (eb < f.ancB), sb1);
                    int fb = md(R[j] - tb1, sb1); int ta = md(f.lbA, f.sa); int fa = f.sa ? md(md(T[j], f.sa) - ta, f.sa) : 0;
                    pred = 1 + comb_cost(fa, f.sa, fb, sb1);
                    fprintf(stderr, "  c2 detail fa=%d sa=%d fb=%d sb1=%d R=%d tb1=%d lbB=%d ancB=%lld eb=%lld T=%d lbA=%d\n", fa, f.sa, fb, sb1, R[j], tb1, f.lbB, f.ancB>>10, eb>>10, T[j], f.lbA); } }
            else if (k == j + 2) { what = "s2"; pred = aug(j, et, T[j], -1, 0, T[j], qi[j].sb ? md(R[j], qi[j].sb) : 0); }
            else { what = "suf"; pred = aug(k - 2, et, T[k - 2], -1, 0); }
            if (pred >= 0 && pred != real) fprintf(stderr, "  k=%d %s real=%d pred=%d\n", k, what, real, pred);
        }
        fprintf(stderr, "  final real=%d old=%d realdelta=%lld bdelta=%d\n", final_rot(s, w), finOld, sumReal + final_rot(s, w) - oldT, b.delta);
    }
};
