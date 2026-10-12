"""Audit of the realizability hierarchy B_bd <= B_walk <= B_real <= B_joint_real <= Delta.

For every interval length L in [4, LMAX] the complete set of 2^L subsets is enumerated, so
Delta(L,n) and Q_real(L,n,u) are exact.  B_bd is the closed form of Theorem 1, B_walk adds
the increment constraint |a_{x+1}-a_x| <= 1 and the endpoint conditions, B_real uses the
exact binary-convolution minimum, and B_joint_real is the LP that combines every window
inequality built from Q_real with the pair budget sum_d r_d = n(n-1)/2.

Outputs
    data/realizability_audit.json          one row per (L,n)
    data/realizability_obstruction.json    the non-tight pairs with certificates
    data/realizability_extended_map.json   exact Delta for L beyond the published map
"""

from __future__ import annotations

import json

import sys

from pathlib import Path

import numpy as np

from common import HERE, save

from joint_bound import B_ph

from realizability_bound import (B_bd, B_from_Q, B_joint_real, B_walk, Q_bd, Q_walk,
                                 enumerate_L, tent_capacity)

def capacity_optimal_profile(L, u, n):
    """A canonical minimiser of sum_x a_x^2 subject to 0 <= a_x <= h_x and sum a_x = n u.

    Boundary windows whose capacity lies below the plateau level k are saturated, the
    plateau holds k, and eta plateau entries are raised to k+1 (bumps placed leftmost).
    """
    h = tent_capacity(L, u)
    M = n * u
    k = 0
    for cand in range(u, -1, -1):
        if cand * (L + u - cand) <= M:
            k = cand
            break
    eta = M - k * (L + u - k)
    a = np.minimum(h, k).astype(np.int64)
    idx = [i for i in range(len(a)) if a[i] == k and h[i] >= k + 1]
    for i in idx[:eta]:
        a[i] = k + 1
    return a, int(k), int(eta), h

def reconstruct_binary(a, L, u):
    """Deterministic reconstruction of the tap vector from a candidate window profile.

    The identity a_x - a_{x-1} = b_{x+u-1} - b_{x-1} determines b uniquely once b_0 is
    fixed, because a_{-(u-1)} = b_0 and, walking left to right,
        b_{x+u-1} = a_x - a_{x-1} + b_{x-1}   whenever both indices exist.
    Windows whose gain index exceeds L-1 carry no new bit and instead impose the
    consistency condition b_{x-1} = a_{x-1} - a_x.  Returns (b, failure) where failure is
    None if the profile is realised by the returned binary vector.
    """
    N = len(a)
    x0 = -(u - 1)
    b = {}
    fail = None
    if int(a[0]) not in (0, 1):
        return None, ("endpoint", x0, "a_{%d} = %d is not 0 or 1, so b_0 does not exist"
                      % (x0, int(a[0])))
    b[0] = int(a[0])
    for i in range(1, N):
        x = x0 + i
        gain = x + u - 1
        loss = x - 1
        if loss >= 0:
            bl = b.get(loss)
            if bl is None:
                return None, ("propagation", x, "b_{%d} is undetermined when required" % loss)
        else:
            bl = 0
        if gain <= L - 1:
            val = int(a[i]) - int(a[i - 1]) + bl
            if val not in (0, 1):
                return None, ("increment", x,
                              "b_{%d} = a_{%d} - a_{%d} + b_{%d} = %d - %d + %d = %d, "
                              "which is not 0 or 1" % (gain, x, x - 1, loss,
                                                       int(a[i]), int(a[i - 1]), bl, val))
            if gain in b and b[gain] != val:
                return None, ("coupling", x,
                              "b_{%d} was fixed to %d earlier but the window at x = %d "
                              "requires %d" % (gain, b[gain], x, val))
            b[gain] = val
        else:
            need = int(a[i - 1]) - int(a[i])
            got = b.get(loss)
            if got is None:
                return None, ("propagation", x, "b_{%d} is undetermined when required" % loss)
            if got != need:
                return None, ("coupling", x,
                              "the window at x = %d has no incoming bit and therefore "
                              "requires b_{%d} = a_{%d} - a_{%d} = %d - %d = %d, whereas the "
                              "left-to-right propagation fixed b_{%d} = %d"
                              % (x, loss, x - 1, x, int(a[i - 1]), int(a[i]), need, loss, got))
    vec = [b.get(i, 0) for i in range(L)]
    if sum(vec) != sum(1 for i in range(L) if b.get(i, 0) == 1):
        pass
    return vec, None

def classify_violation(a, L, u):
    """First constraint that the capacity-optimal profile violates."""
    h = tent_capacity(L, u)
    N = len(a)
    if int(a[0]) > 1:
        return ("endpoint", 0,
                "a_{-(u-1)} = %d but a_{-(u-1)} = b_0 is 0 or 1" % int(a[0]))
    if int(a[N - 1]) > 1:
        return ("endpoint", N - 1,
                "a_{L-1} = %d but a_{L-1} = b_{L-1} is 0 or 1" % int(a[N - 1]))
    for i in range(N - 1):
        d = int(a[i + 1]) - int(a[i])
        if abs(d) >= 2:
            return ("increment", i,
                    "a_{x+1} - a_x = %d at x = %d, but a_{x+1} - a_x = b_{x+u} - b_x lies in {-1,0,1}"
                    % (d, i - (u - 1)))
    return ("coupling", None,
            "increments and endpoints are admissible, yet no binary vector reproduces the "
            "profile: the shared-bit coupling between the down-step at x and the up-step "
            "at x+u fails")

def profile_certificate(a, L, u):
    """Human-checkable certificate: deterministic reconstruction plus increment check."""
    vec, fail = reconstruct_binary(a, L, u)
    kind, pos, why = classify_violation(a, L, u)
    if fail is None:
        return dict(realizable=True, kind="-", position=None,
                    detail="realised by b = (%s)" % ",".join(str(v) for v in vec),
                    witness="".join(str(v) for v in vec))
    return dict(realizable=False, kind=fail[0], position=fail[1], detail=fail[2],
                coarse_kind=kind, witness=None)

def B_joint_int(L, n, Qfun, us=None, time_limit=20):
    """min m over INTEGER r_d in [0,m] with sum r_d = C and every window inequality."""
    from scipy.optimize import milp, LinearConstraint, Bounds
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
    A = np.zeros((len(rows) + (L - 1), L))
    lb = np.zeros(len(rows) + (L - 1))
    ub = np.zeros(len(rows) + (L - 1))
    for i, (coef, b) in enumerate(zip(rows, rhs)):
        A[i, :] = coef
        lb[i] = b
        ub[i] = np.inf
    for d in range(1, L):
        A[len(rows) + d - 1, d] = 1.0
        A[len(rows) + d - 1, 0] = -1.0
        lb[len(rows) + d - 1] = 0.0
        ub[len(rows) + d - 1] = 0.0
    A_eq = np.zeros((1, L))
    A_eq[0, 1:L] = 1.0
    c = np.zeros(L)
    c[0] = 1.0
    res = milp(c=c, constraints=[LinearConstraint(A, lb, ub),
                                 LinearConstraint(A_eq, float(C), float(C))],
               integrality=np.ones(L), bounds=Bounds(np.zeros(L), np.full(L, float(C))),
               options=dict(time_limit=time_limit))
    if not res.success or res.x is None:
        return None
    return int(np.ceil(res.x[0] - 1e-9))

def joint_int_feasible(L, n, m, Qfun, us=None, time_limit=60):
    """Is there an integer vector r with 0 <= r_d <= m, sum r_d = C, satisfying every
    window inequality built from Qfun?  Returns (feasible, r) or (False, None)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    C = n * (n - 1) // 2
    us = list(range(2, L + 1)) if us is None else us
    rows, rhs = [], []
    for u in us:
        q = Qfun(L, u, n)
        if q is None:
            continue
        need = q - n * u
        if need <= 0:
            continue
        coef = np.zeros(L - 1)
        for d in range(1, u):
            coef[d - 1] = 2.0 * (u - d)
        rows.append(coef)
        rhs.append(float(need))
    if not rows:
        rows = [np.zeros(L - 1)]
        rhs = [0.0]
    A = np.array(rows)
    c = np.zeros(L - 1)
    res = milp(c=c, constraints=[LinearConstraint(A, rhs, np.inf),
                                 LinearConstraint(np.ones((1, L - 1)), float(C), float(C))],
               integrality=np.ones(L - 1),
               bounds=Bounds(np.zeros(L - 1), np.full(L - 1, float(min(m, C)))),
               options=dict(time_limit=time_limit))
    if not res.success or res.x is None:
        return None, None
    r = np.round(res.x).astype(int)
    ok = (r.sum() == C) and (r.min() >= 0) and (r.max() <= m)
    if A is not None and ok:
        ok = bool((A @ r >= np.array(rhs) - 1e-6).all())
    return (bool(ok), r.tolist()) if ok else (False, None)

def per_u_table(L, n, Q):
    out = []
    for u in range(2, L + 1):
        qb = Q_bd(L, u, n)
        qr = Q[(n, u)]
        cb = -((-(qb - n * u)) // (u * (u - 1))) if qb - n * u > 0 else 0
        cr = -((-(qr - n * u)) // (u * (u - 1))) if qr - n * u > 0 else 0
        out.append(dict(u=u, Q_bd=qb, Q_real=qr, improvement=qr - qb,
                        bound_from_Q_bd=cb, bound_from_Q_real=cr))
    return out

def finite_b2g(all_delta, gmax=8):
    """F_g(L) = largest n with Delta(L,n) <= g = largest size of a positive-difference-bounded set in L points."""
    out = []
    for L in sorted(all_delta):
        D = all_delta[L]
        row = dict(L=L)
        for g in range(1, gmax + 1):
            cand = [n for n in D if D[n] < 10 ** 8 and D[n] <= g]
            row[f"F_{g}"] = max(cand) if cand else 0
        out.append(row)
    return out
