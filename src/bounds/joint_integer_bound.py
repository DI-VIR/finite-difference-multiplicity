"""Integer joint window-pair relaxation.

Variables are the difference multiplicities r_d, d = 1..L-1, and the maximum m.
Constraints, all valid for every tap set of cardinality n in I_L:

   (C1) sum_d r_d = n(n-1)/2                                   (every unordered pair)
   (C2) 0 <= r_d <= m                                           (definition of m)
   (C3) 2 sum_{d<u} (u-d) r_d >= Q_bd(u) - n u   for every u    (window identity
                                                                 plus the boundary
                                                                 capacity minimum)

The integer optimum m* is a valid lower bound on Delta(L,n).  It dominates B_bd and
B_ph, both of which are the two special cases obtained by keeping only (C2)+(C3) or only
(C1)+(C2) respectively. On tested pairs it equals max(B_bd,B_ph);
it does not certify the remaining exact-map gaps.
"""

from __future__ import annotations

import json

import sys

import numpy as np

from scipy.optimize import milp, LinearConstraint, Bounds

from cocapacity_bound import capacity_min, bound as bd_bound

def Q_bd_u(L, u, M):
    return capacity_min(L, u, M, u)

def joint_integer(L, n, m_max=None, time_limit=60):
    """Smallest integer m for which (C1)-(C3) are feasible; None if infeasible up to m_max."""
    C = n * (n - 1) // 2
    m_max = m_max or n - 1
    rows, rhs = [], []
    for u in range(2, L + 1):
        need = Q_bd_u(L, u, n * u) - n * u
        if need <= 0:
            continue
        coef = [0.0] * L
        for d in range(1, u):
            coef[d] -= 2.0 * (u - d)
        rows.append(coef)
        rhs.append(-float(need))
    A = np.array(rows, dtype=float) if rows else np.zeros((0, L))
    b = np.array(rhs, dtype=float) if rows else np.zeros(0)
    for m in range(1, m_max + 1):
        A_ub = np.vstack([A] + [np.eye(L)[d] - np.eye(L)[0] for d in range(1, L)]) \
            if rows else np.array([np.eye(L)[d] - np.eye(L)[0] for d in range(1, L)])
        b_ub = np.concatenate([b, np.zeros(L - 1)])
        c = np.zeros(L)
        integrality = np.ones(L)
        bounds = Bounds(np.array([float(m)] + [0.0] * (L - 1)),
                        np.array([float(m)] + [float(min(m, n - 1))] * (L - 1)))
        res = milp(c=c, constraints=[LinearConstraint(A_ub, -np.inf, b_ub),
                                     LinearConstraint(np.eye(L)[1:].sum(0).reshape(1, L), C, C)],
                   integrality=integrality, bounds=bounds,
                   options=dict(time_limit=time_limit, mip_rel_gap=0.0))
        if res.success:
            r=np.rint(res.x[1:]).astype(int)
            assert int(r.sum())==C and all(0<=v<=m for v in r)
            assert all(2*sum((u-d)*int(r[d-1]) for d in range(1,u))>=Q_bd_u(L,u,n*u)-n*u
                       for u in range(2,L+1))
            return m,r
        if res.status!=2:
            return None,None
    return None, None
