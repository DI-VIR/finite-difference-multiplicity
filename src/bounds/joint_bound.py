"""Joint window-pair relaxation: a lower bound that uses the window identity and the
exact pair budget simultaneously.

For every tap set the difference multiplicities r_d = r_T(d) satisfy two families of
facts that the boundary-capacity bound uses only separately:

  (i)  window identity, for every u:  2 * sum_{d<u} (u-d) r_d  =  sum_x a_x(u)^2 - n u
       and sum_x a_x(u)^2 >= Q_bd(u), hence
                    2 * sum_{d<u} (u-d) r_d >= Q_bd(u) - n u;

  (ii) pair budget:  sum_{d=1}^{L-1} r_d = n(n-1)/2.

Minimising max_d r_d subject to both gives a valid lower bound on Delta(L,n) that is
never weaker than either B_bd or B_ph:

       B_joint(L,n) = ceil( min  m   s.t.  0 <= r_d <= m,  sum_d r_d = n(n-1)/2,
                                          2 sum_{d<u}(u-d) r_d >= Q_bd(u) - n u  for all u )

The minimum over real r is a relaxation of the integer problem, so its ceiling is still
a valid integer lower bound.  B_joint >= B_bd and B_joint >= B_ph both follow because
the constraints alone imply max_d r_d >= (Q_bd(u)-nu)/(u(u-1)) for every u and
max_d r_d >= n(n-1)/(2(L-1)).
"""

from __future__ import annotations

import json

import sys

import numpy as np

from scipy.optimize import linprog

def capacity_min_square_sum(L, u, M):
    """Q_bd(u) of Theorem 1."""
    k = 0
    for cand in range(u, -1, -1):
        if cand * (L + u - cand) <= M:
            k = cand
            break
    eta = M - k * (L + u - k)
    return ((k - 1) * k * (2 * k - 1)) // 3 + (L + u - 2 * k + 1) * k * k + (2 * k + 1) * eta

def B_bd(L, n):
    best = 0
    for u in range(2, L + 1):
        M = n * u
        D = u * (u - 1)
        val = capacity_min_square_sum(L, u, M) - M
        if val > 0:
            ceil = -((-val) // D)
            if ceil > best:
                best = ceil
    return best

def B_ph(L, n):
    val = n * (n - 1) // 2
    return -((-val) // (L - 1))

def B_joint(L, n, us=None, return_raw=False):
    """LP relaxation of the joint window-pair problem; returns ceil of the LP value."""
    C = n * (n - 1) // 2
    if C == 0:
        return 0, 0.0
    us = list(range(2, L + 1)) if us is None else us
    rows = []
    rhs = []
    for u in us:
        need = capacity_min_square_sum(L, u, n * u) - n * u
        if need <= 0:
            continue
        coef = np.zeros(L)
        for d in range(1, u):
            coef[d] -= 2.0 * (u - d)
        rows.append(coef)
        rhs.append(-float(need))
    if not rows:
        return B_ph(L, n), float(B_ph(L, n))
    # variables: r_1..r_{L-1} (index 0 unused), then m
    A_ub = np.zeros((len(rows) + (L - 1), L))
    b_ub = np.zeros(len(rows) + (L - 1))
    for i, (coef, b) in enumerate(zip(rows, rhs)):
        A_ub[i, :] = coef
        b_ub[i] = b
    for d in range(1, L):
        A_ub[len(rows) + d - 1, d] = 1.0
        A_ub[len(rows) + d - 1, 0] = -1.0   # r_d <= m  (m stored in slot 0)
        b_ub[len(rows) + d - 1] = 0.0
    A_eq = np.zeros((1, L))
    A_eq[0, 1:L] = 1.0
    b_eq = np.array([float(C)])
    c = np.zeros(L)
    c[0] = 1.0
    bounds = [(0, None)] * L
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    if not res.success:
        return None, None
    raw = float(res.x[0])
    value = int(np.ceil(raw - 1e-9))
    return (value, raw) if return_raw else (value, raw)

def table(pairs):
    out = []
    for L, n in pairs:
        bb = B_bd(L, n)
        bp = B_ph(L, n)
        bj, raw = B_joint(L, n)
        out.append(dict(L=L, n=n, B_ph=bp, B_bd=bb, B_joint=bj, lp_value=round(raw, 4),
                        improves_bd=bj - bb, improves_max=max(0, bj - max(bb, bp))))
    return out
