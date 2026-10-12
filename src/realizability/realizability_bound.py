"""Realizability hierarchy for the finite-interval overlap problem.

Notation.  T is an n-subset of {0,...,L-1}; r_d = #{(a,b) in T^2 : b-a = d} for d >= 1;
Delta(L,n) = min over T of max_{d>=1} r_d.  For a window length u put N = L+u-1, let
window x run over x = -(u-1),...,L-1 and let

        a_x = #{t in T : x <= t <= x+u-1},      h_x = #{j in [0,L-1] : x <= j <= x+u-1}.

Then sum_x a_x = n u and the window identity holds:

        sum_x a_x^2 = n u + 2 * sum_{d=1}^{u-1} (u-d) r_d.                        (W)

Every relaxation of the set of attainable profiles a = (a_x) therefore produces a lower
bound on Delta: if Q(u) <= min sum_x a_x^2 over the relaxed set, then for any T,

        max_d r_d  >=  ( Q(u) - n u ) / ( u (u-1) ),
        B(L,n)     =  max_{u>=2} ceil( ( Q(u) - n u ) / ( u (u-1) ) ).

Three nested relaxations are considered here.

  Q_bd    : 0 <= a_x <= h_x, a integer, sum a_x = n u.  (Theorem 1 of the paper, closed form.)
  Q_walk  : Q_bd plus the increment constraint |a_{x+1} - a_x| <= 1 together with the
            endpoint conditions a_{-(u-1)} = b_0 in {0,1} and a_{L-1} = b_{L-1} in {0,1}.
            Both are implied by a_{x+1} - a_x = b_{x+u} - b_x with b binary, so Q_walk is
            still a relaxation.  Computed by dynamic programming.
  Q_real  : the exact minimum over binary b with popcount n.  Computed by complete
            enumeration of all 2^L subsets for L <= LMAX, and by exact MILP otherwise.

Since the feasible sets are nested, Q_bd <= Q_walk <= Q_real, hence
B_bd <= B_walk <= B_real <= Delta.  B_joint_real additionally imposes the pair budget
sum_d r_d = n(n-1)/2 together with every window inequality built from Q_real.
"""

from __future__ import annotations

import json

import sys

from pathlib import Path

import numpy as np

from common import HERE, save

from joint_bound import B_ph, capacity_min_square_sum

def tent_capacity(L, u):
    """h_x for x = -(u-1),...,L-1 (length N = L+u-1)."""
    x = np.arange(-(u - 1), L, dtype=np.int64)
    lo = np.maximum(x, 0)
    hi = np.minimum(x + u - 1, L - 1)
    return np.maximum(hi - lo + 1, 0)

def Q_bd(L, u, n):
    return capacity_min_square_sum(L, u, n * u)

def Q_walk(L, u, n):
    """Minimum of sum_x a_x^2 over integer profiles with cumulative sum n u, capacity
    constraints, |increment| <= 1 and endpoints in {0,1}.  None if infeasible."""
    h = tent_capacity(L, u)
    N = len(h)
    M = n * u
    INF = float("inf")
    # dp[s, a] : minimum cost after the current window, total s, current value a
    dp = np.full((M + 1, u + 1), INF)
    for a0 in (0, 1):
        if a0 <= h[0] and a0 <= M:
            dp[a0, a0] = float(a0 * a0)
    for i in range(1, N):
        cap = int(min(h[i], u))
        ndp = np.full((M + 1, u + 1), INF)
        for a in range(cap + 1):
            best = None
            for da in (-1, 0, 1):
                p = a + da
                if p < 0 or p > u:
                    continue
                col = dp[:, p]
                if a == 0:
                    cand = col
                else:
                    cand = np.full(M + 1, INF)
                    if M + 1 - a > 0:
                        cand[a:] = col[: M + 1 - a]
                best = cand if best is None else np.minimum(best, cand)
            if best is None:
                continue
            ndp[:, a] = best + a * a
        if i == N - 1:                      # last window: a_{L-1} = b_{L-1} in {0,1}
            ndp[:, 2:] = INF
        dp = ndp
        if not np.isfinite(dp[:, : min(2, u + 1)]).any() and i < N - 1:
            pass
    val = dp[M, 0] if u == 0 else min(dp[M, 0], dp[M, 1] if u >= 1 else INF)
    return int(val) if np.isfinite(val) else None

def enumerate_L(L, chunk=1 << 18, with_evidence=False):
    """Complete enumeration of all 2^L subsets of [0,L-1].

    Returns (Delta_by_n, Qreal_by_n_u) with
        Delta_by_n[n]      = min over subsets of size n of max_d r_d   (exact)
        Qreal_by_n_u[n][u] = n u + 2 * min over subsets of size n of sum_{d<u}(u-d) r_d
    Both are exact; no sampling and no relaxation is involved.
    """
    total = 1 << L
    delta_min = np.full(L + 1, 10 ** 9, dtype=np.int64)
    w_min = np.full((L + 1, L + 1), 10 ** 9, dtype=np.int64)   # [n, u] -> min sum_{d<u}(u-d) r_d
    optimum_counts = [0]*(L+1)
    anchored_counts = [0]*(L+1)
    class_parts = [[] for _ in range(L+1)]
    witnesses = [None]*(L+1)
    witness_keys = [-1]*(L+1)
    start = 0
    while start < total:
        end = min(start + chunk, total)
        s = np.arange(start, end, dtype=np.uint64)
        cnt = np.bitwise_count(s).astype(np.int64)
        R = np.empty((end - start, L), dtype=np.int64)
        for d in range(1, L):
            R[:, d - 1] = np.bitwise_count(s & (s >> np.uint64(d)))
        R[:, L - 1] = 0                      # r_L does not exist; keeps the shape uniform
        # W_u = sum_{d<u} (u-d) r_d  equals the double cumulative sum of (r_1,...,r_{L-1})
        B = np.cumsum(np.cumsum(R, axis=1), axis=1)            # B[:, u-2] = W_u for u >= 2
        W = np.zeros((end - start, L + 1), dtype=np.int64)
        W[:, 2:] = B[:, : L - 1]
        dmax = R.max(axis=1)
        order = np.argsort(cnt, kind="stable")
        cs = cnt[order]
        starts = np.searchsorted(cs, np.arange(L + 1))
        ends = np.searchsorted(cs, np.arange(L + 1), side="right")
        for n in range(L + 1):
            if starts[n] >= ends[n]:
                continue
            sl = order[starts[n]:ends[n]]
            current = int(dmax[sl].min())
            if with_evidence and 2 <= n <= L and current <= delta_min[n]:
                if current < delta_min[n]:
                    optimum_counts[n] = anchored_counts[n] = 0
                    class_parts[n] = []
                    witnesses[n] = None
                    witness_keys[n] = -1
                optimal = s[sl[dmax[sl] == current]]
                optimum_counts[n] += len(optimal)
                anchored = optimal[(optimal & 1) != 0]
                anchored_counts[n] += len(anchored)
                if len(anchored):
                    reverse = np.zeros(len(anchored), dtype=np.uint64)
                    for bit in range(L):
                        reverse = (reverse << np.uint64(1)) | ((anchored >> np.uint64(bit)) & np.uint64(1))
                    span = np.floor(np.log2(anchored)).astype(np.uint64)
                    reflected = reverse >> (np.uint64(L-1)-span)
                    class_parts[n].append(np.minimum(anchored,reflected))
                    i = int(reverse.argmax())
                    if int(reverse[i]) > witness_keys[n]:
                        witness_keys[n] = int(reverse[i])
                        witnesses[n] = [bit for bit in range(L) if (int(anchored[i]) >> bit) & 1]
            delta_min[n] = min(delta_min[n], current)
            w_min[n] = np.minimum(w_min[n], W[sl].min(axis=0))
        start = end
    qreal = {}
    for n in range(L + 1):
        for u in range(2, L + 1):
            qreal[(n, u)] = int(n * u + 2 * w_min[n, u])
    values = {n: int(delta_min[n]) for n in range(L + 1)}
    if with_evidence:
        evidence = {n: dict(witness=witnesses[n], all_optimal_subsets=optimum_counts[n],
                           anchored_optima=anchored_counts[n],
                           translation_reflection_classes=len(np.unique(np.concatenate(class_parts[n])))
                           if class_parts[n] else 0) for n in range(2,L+1)}
        return values,qreal,evidence
    return values,qreal

def B_from_Q(L, n, Qfun, us=None):
    best, arg = 0, None
    for u in range(2, L + 1) if us is None else us:
        q = Qfun(L, u, n)
        if q is None:
            continue
        val = q - n * u
        if val <= 0:
            continue
        c = -((-val) // (u * (u - 1)))
        if c > best:
            best, arg = c, u
    return best, arg

def B_bd(L, n):
    return B_from_Q(L, n, Q_bd)

def B_walk(L, n, cache=None, us=None):
    def f(L, u, n):
        key = (L, u, n)
        if cache is not None and key in cache:
            return cache[key]
        v = Q_walk(L, u, n)
        if cache is not None:
            cache[key] = v
        return v
    return B_from_Q(L, n, f, us=us)

def B_joint_real(L, n, Qfun, us=None):
    """min m s.t. 0 <= r_d <= m, sum_d r_d = C, 2 sum_{d<u}(u-d) r_d >= Q(u) - n u."""
    from scipy.optimize import linprog
    C = n * (n - 1) // 2
    if C == 0:
        return 0
    us = list(range(2, L + 1)) if us is None else us
    rows, rhs = [], []
    for u in us:
        q = Qfun(L, u, n)
        if q is None:
            continue
        need = q - n * u
        if need <= 0:
            continue
        coef = np.zeros(L)
        for d in range(1, u):
            coef[d] -= 2.0 * (u - d)
        rows.append(coef)
        rhs.append(-float(need))
    if not rows:
        return B_ph(L, n)
    A_ub = np.zeros((len(rows) + (L - 1), L))
    b_ub = np.zeros(len(rows) + (L - 1))
    for i, (coef, b) in enumerate(zip(rows, rhs)):
        A_ub[i, :] = coef
        b_ub[i] = b
    for d in range(1, L):
        A_ub[len(rows) + d - 1, d] = 1.0
        A_ub[len(rows) + d - 1, 0] = -1.0
    A_eq = np.zeros((1, L))
    A_eq[0, 1:L] = 1.0
    c = np.zeros(L)
    c[0] = 1.0
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=np.array([float(C)]),
                  bounds=[(0, None)] * L, method="highs")
    if not res.success:
        return None
    return int(np.ceil(res.x[0] - 1e-9))
