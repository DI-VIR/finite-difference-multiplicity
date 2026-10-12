"""Verify the scope of the joint window-pair relaxation on every enumerated pair.

For each instance we compare

    B_bd                        (Theorem 1, window identity at a single u)
    B_ph                        (pair counting)
    B_joint^Z                   (Theorem 4, both together, integer feasible set)
    B_bd2                       (complement-capacity variant, valid but not adopted)

and record whether the joint relaxation attains max(B_bd, B_ph), i.e. whether an integer
difference-multiplicity vector satisfying every window constraint exists at that level.
"""

from __future__ import annotations

import csv

import json

import sys

import numpy as np

from scipy.optimize import milp, LinearConstraint, Bounds

from common import HERE, save

from cocapacity_bound import capacity_min, bound as bd_bound

def Q_bd_u(L, u, M):
    return capacity_min(L, u, M, u)

def B_ph(L, n):
    return -((-n * (n - 1) // 2) // (L - 1))

def feasible_at(L, n, m, time_limit=20):
    C = n * (n - 1) // 2
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
    ub_rows = [np.eye(L)[d] - np.eye(L)[0] for d in range(1, L)]
    A_ub = np.array(rows + ub_rows, dtype=float)
    b_ub = np.array(rhs + [0.0] * (L - 1), dtype=float)
    eye = np.zeros((1, L))
    eye[0, 1:] = 1.0
    res = milp(c=np.zeros(L),
               constraints=[LinearConstraint(A_ub, -np.inf, b_ub), LinearConstraint(eye, C, C)],
               integrality=np.ones(L),
               bounds=Bounds(np.array([float(m)] + [0.0] * (L - 1)),
                             np.array([float(m)] + [float(min(m, n - 1))] * (L - 1))),
               options=dict(time_limit=time_limit, mip_rel_gap=0.0))
    if res.success:
        r=np.rint(res.x[1:]).astype(int).tolist()
        assert all(0<=v<=m for v in r) and sum(r)==C
        assert all(2*sum((u-d)*r[d-1] for d in range(1,u))>=Q_bd_u(L,u,n*u)-n*u
                   for u in range(2,L+1))
        return "feasible",r
    if res.status == 2:
        return "infeasible_proved",None
    return "unknown",None
