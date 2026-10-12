"""Exact weighted-pair model for one binary-window square-sum minimum.

Only the small validation instances below have been executed in this revision.
The formulation carries no large-instance runtime guarantee.
"""
from itertools import combinations
from pathlib import Path
import json
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]
VALIDATION_INSTANCES = (
    (4, 2, 2), (4, 2, 4), (6, 3, 2), (6, 3, 6),
    (8, 4, 3), (8, 4, 8), (10, 4, 5), (10, 6, 10),
    (11, 5, 6), (12, 5, 2), (12, 6, 7), (12, 6, 12),
)


def square_sum(marks, u):
    return len(marks) * u + 2 * sum(max(0, u - (b - a))
                                  for a, b in combinations(sorted(marks), 2))


def optimize(L, n, u, time_limit=60.0):
    if not (isinstance(L, int) and isinstance(n, int) and isinstance(u, int)
            and 2 <= n <= L and 2 <= u <= L):
        raise ValueError('Require integer 2 <= n,u <= L.')
    if time_limit <= 0:
        raise ValueError('time_limit must be positive.')
    products = [(i, d) for d in range(1, u) for i in range(L - d)]
    size = L + len(products)
    objective = np.zeros(size)
    row, col, val = [], [], []
    lower, upper = [float(n)], [float(n)]
    for i in range(L):
        row.append(0); col.append(i); val.append(1.0)
    for j, (i, d) in enumerate(products, L):
        objective[j] = u - d
        for entries, lo, hi in (
            (((j, 1), (i, -1)), -np.inf, 0),
            (((j, 1), (i + d, -1)), -np.inf, 0),
            (((j, 1), (i, -1), (i + d, -1)), -1, np.inf),
        ):
            r = len(lower)
            for c, v in entries:
                row.append(r); col.append(c); val.append(float(v))
            lower.append(lo); upper.append(hi)
    matrix = coo_matrix((val, (row, col)), shape=(len(lower), size)).tocsc()
    started = time.perf_counter()
    result = milp(objective, integrality=np.r_[np.ones(L), np.zeros(len(products))],
                  bounds=Bounds(np.zeros(size), np.ones(size)),
                  constraints=LinearConstraint(matrix, lower, upper),
                  options={'time_limit': float(time_limit), 'mip_rel_gap': 0.0})
    seconds = round(time.perf_counter() - started, 6)
    if result.status != 0 or result.x is None:
        raise RuntimeError('No certified solver optimum; status=%s: %s'
                           % (result.status, result.message))
    bits = np.rint(result.x[:L]).astype(int)
    if not np.allclose(result.x[:L], bits, atol=1e-7, rtol=0):
        raise ArithmeticError('Nonintegral endpoint solution.')
    marks = [int(i) for i in np.flatnonzero(bits)]
    if len(marks) != n or not all(v in (0, 1) for v in bits):
        raise ArithmeticError('Invalid binary witness.')
    q = square_sum(marks, u)
    if abs(result.fun - (q - n * u) / 2) > 1e-6:
        raise ArithmeticError('Solver objective disagrees with exact pair counting.')
    return dict(L=L, n=n, u=u, Q_real=q, witness=marks, solver_status=int(result.status),
                seconds=seconds, mip_gap=float(result.mip_gap),
                mip_nodes=int(result.mip_node_count), time_limit=float(time_limit))


def validate():
    records = []
    for L, n, u in VALIDATION_INSTANCES:
        record = optimize(L, n, u)
        exhaustive = min(square_sum(t, u) for t in combinations(range(L), n))
        assert record['Q_real'] == exhaustive, (L, n, u)
        record['independent_enumeration_Q'] = exhaustive
        record['enumerated_sets'] = sum(1 for _ in combinations(range(L), n))
        records.append(record)
    return dict(status='PASS', instances=len(records), records=records,
                comparison='Independent unordered-pair enumeration',
                large_instance_runs=False, scalability_claim=False)


if __name__ == '__main__':
    report = validate()
    (ROOT / 'data/binary_window_milp_validation.json').write_text(
        json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
