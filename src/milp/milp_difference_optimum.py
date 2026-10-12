from __future__ import annotations

import argparse

import csv

import json

import math

import time

from collections import Counter

from pathlib import Path

from typing import Iterable

import numpy as np

from scipy.optimize import Bounds, LinearConstraint, milp

from scipy.sparse import coo_matrix

DEFAULT_CASES = [
    (24, 11, "claimed_missing"),
    (32, 15, "claimed_existing_s4"),
    (40, 19, "claimed_missing"),
    (48, 22, "claimed_missing"),
    (80, 14, "claimed_missing"),
    (64, 29, "extended"),
    (80, 37, "extended"),
]

def integer_cs_bound(length: int, n_taps: int) -> tuple[int, int]:
    best = 0
    best_u = 2
    for u in range(2, length + 1):
        window_count = length + u - 1
        total = n_taps * u
        q, s = divmod(total, window_count)
        min_square_sum = (window_count - s) * q * q + s * (q + 1) * (q + 1)
        bound = math.ceil((min_square_sum - total) / (u * (u - 1)))
        if bound > best:
            best = bound
            best_u = u
    return best, best_u

def parse_cases(text: str) -> list[tuple[int, int, str]]:
    cases: list[tuple[int, int, str]] = []
    for raw in text.split(","):
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split(":")
        if len(parts) == 2:
            length, n_taps = map(int, parts)
            tag = "custom"
        elif len(parts) == 3:
            length, n_taps = map(int, parts[:2])
            tag = parts[2]
        else:
            raise ValueError(f"bad case item: {raw}")
        cases.append((length, n_taps, tag))
    return cases

def pair_index(length: int) -> tuple[dict[tuple[int, int], int], dict[int, list[int]], int]:
    offset = length
    pair_to_var: dict[tuple[int, int], int] = {}
    diff_to_vars: dict[int, list[int]] = {d: [] for d in range(1, length)}
    cursor = offset
    for i in range(length):
        for j in range(i + 1, length):
            pair_to_var[(i, j)] = cursor
            diff_to_vars[j - i].append(cursor)
            cursor += 1
    return pair_to_var, diff_to_vars, cursor

def build_constraints(length: int, n_taps: int, max_delta: int) -> tuple[LinearConstraint, Bounds, np.ndarray]:
    pair_to_var, diff_to_vars, n_vars = pair_index(length)
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    lower: list[float] = []
    upper: list[float] = []

    def add_row(items: Iterable[tuple[int, float]], lb: float, ub: float) -> None:
        row = len(lower)
        for col, val in items:
            rows.append(row)
            cols.append(col)
            vals.append(val)
        lower.append(lb)
        upper.append(ub)

    # Select exactly n tap positions.
    add_row(((i, 1.0) for i in range(length)), float(n_taps), float(n_taps))

    # Translation symmetry break. Every feasible difference profile can be
    # shifted so that the smallest selected tap is at position 0.
    add_row([(0, 1.0)], 1.0, 1.0)

    # Pair variables y_ij represent x_i AND x_j.
    for (i, j), y in pair_to_var.items():
        add_row([(y, 1.0), (i, -1.0)], -np.inf, 0.0)
        add_row([(y, 1.0), (j, -1.0)], -np.inf, 0.0)
        add_row([(y, -1.0), (i, 1.0), (j, 1.0)], -np.inf, 1.0)

    # Difference multiplicity upper bounds.
    for d in range(1, length):
        add_row(((var, 1.0) for var in diff_to_vars[d]), -np.inf, float(max_delta))

    matrix = coo_matrix((vals, (rows, cols)), shape=(len(lower), n_vars)).tocsr()
    constraints = LinearConstraint(matrix, np.array(lower), np.array(upper))
    bounds = Bounds(np.zeros(n_vars), np.ones(n_vars))
    integrality = np.ones(n_vars, dtype=np.int8)
    return constraints, bounds, integrality

def counts_from_taps(taps: list[int], length: int) -> list[int]:
    tap_set = set(taps)
    return [0] + [sum(1 for t in taps if t + d in tap_set) for d in range(1, length)]

def verify_taps(taps: list[int], length: int, n_taps: int, max_delta: int) -> tuple[bool, int, list[int], dict[int, int]]:
    if len(taps) != n_taps or len(set(taps)) != n_taps:
        return False, -1, [], {}
    if not taps or min(taps) < 0 or max(taps) >= length:
        return False, -1, [], {}
    counts = counts_from_taps(taps, length)
    delta = max(counts[1:]) if length > 1 else 0
    peaks = [d for d in range(1, length) if counts[d] == delta]
    histogram = dict(sorted(Counter(counts[1:]).items(), reverse=True))
    return delta <= max_delta, delta, peaks, histogram

def solve_decision(length: int, n_taps: int, max_delta: int, time_limit: float) -> dict[str, object]:
    start = time.perf_counter()
    constraints, bounds, integrality = build_constraints(length, n_taps, max_delta)
    c = np.zeros(len(integrality), dtype=float)
    result = milp(
        c=c,
        integrality=integrality,
        bounds=bounds,
        constraints=constraints,
        options={"time_limit": time_limit, "mip_rel_gap": 0.0, "disp": False},
    )
    elapsed = time.perf_counter() - start
    status = int(result.status)
    message = str(result.message)
    feasible = bool(result.success and result.x is not None)
    taps: list[int] = []
    verified = False
    delta = ""
    peaks: list[int] = []
    histogram: dict[int, int] = {}
    if feasible:
        x = np.asarray(result.x[:length])
        taps = [i for i, value in enumerate(x) if value >= 0.5]
        verified, found_delta, peaks, histogram = verify_taps(taps, length, n_taps, max_delta)
        delta = found_delta
        feasible = feasible and verified
    node_count = getattr(result, "mip_node_count", "")
    mip_gap = getattr(result, "mip_gap", "")
    dual_bound = getattr(result, "mip_dual_bound", "")
    return {
        "target_delta": max_delta,
        "scipy_status": status,
        "solver_success": bool(result.success),
        "solver_message": message,
        "decision": "feasible" if feasible else ("infeasible" if status == 2 else "undetermined"),
        "verified_solution": verified,
        "found_delta": delta,
        "taps": taps,
        "peak_d": peaks,
        "count_histogram": histogram,
        "elapsed_seconds": elapsed,
        "mip_node_count": node_count,
        "mip_gap": mip_gap,
        "mip_dual_bound": dual_bound,
        "variables": len(integrality),
        "constraints": int(constraints.A.shape[0]) if hasattr(constraints, "A") else "",
    }

def certify_case(
    length: int,
    n_taps: int,
    tag: str,
    time_limit_per_target: float,
    max_extra_delta: int,
) -> dict[str, object]:
    lower, best_u = integer_cs_bound(length, n_taps)
    attempts = []
    optimum = ""
    taps: list[int] = []
    peaks: list[int] = []
    histogram: dict[int, int] = {}
    case_start = time.perf_counter()
    for target in range(lower, min(n_taps, lower + max_extra_delta + 1)):
        attempt = solve_decision(length, n_taps, target, time_limit_per_target)
        attempts.append(attempt)
        if attempt["decision"] == "feasible":
            optimum = target
            taps = list(attempt["taps"])
            peaks = list(attempt["peak_d"])
            histogram = dict(attempt["count_histogram"])
            break
        if attempt["decision"] == "undetermined":
            break
    completed = optimum != "" and all(a["decision"] == "infeasible" for a in attempts[:-1])
    return {
        "L": length,
        "n": n_taps,
        "tag": tag,
        "B_int": lower,
        "u_star": best_u,
        "status": "optimal_certified" if completed else "not_certified",
        "delta_star": optimum,
        "certification_logic": (
            "MILP infeasibility for all tested targets below delta_star, and MILP feasible verified solution at delta_star."
            if completed and attempts[:-1]
            else "Exact integer-CS lower bound equals the first MILP-feasible target."
            if completed
            else "Solver did not complete the required infeasibility/feasibility chain."
        ),
        "attempt_count": len(attempts),
        "attempts": attempts,
        "taps": taps,
        "peak_d": peaks,
        "count_histogram": histogram,
        "elapsed_seconds": time.perf_counter() - case_start,
    }

def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
