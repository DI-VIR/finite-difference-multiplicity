"""Finite-interval bounds and exhaustive small tap-set optima.

All bounds are integer expressions; no floating-point ceiling is used.
"""

from __future__ import annotations

from collections import defaultdict

from pathlib import Path

import argparse

import csv

import json

import math

import time

def ceil_div(a: int, b: int) -> int:
    return -(-a // b)

def min_boundary_square(length: int, n_taps: int, u: int) -> tuple[int, int, int]:
    total = n_taps * u
    max_layer = min(n_taps, u)
    lo, hi = 0, max_layer
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid * (length + u - mid) <= total:
            lo = mid
        else:
            hi = mid - 1
    k = lo
    rem = total - k * (length + u - k)
    square = (k - 1) * k * (2 * k - 1) // 3
    square += (length + u - 2 * k + 1) * k * k + (2 * k + 1) * rem
    return square, k, rem

def bounds(length: int, n_taps: int) -> dict[str, int]:
    if not 2 <= n_taps <= length:
        raise ValueError((length, n_taps))
    pig = ceil_div(n_taps * (n_taps - 1), 2 * (length - 1))
    cs = inte = bd = 0
    witness = 2
    for u in range(2, length + 1):
        num_windows = length + u - 1
        total = n_taps * u
        denom = u * (u - 1)
        cs = max(cs, ceil_div(total * total - total * num_windows, num_windows * denom))
        q, rem = divmod(total, num_windows)
        int_square = num_windows * q * q + rem * (2 * q + 1)
        inte = max(inte, ceil_div(int_square - total, denom))
        bd_square, _, _ = min_boundary_square(length, n_taps, u)
        trial = ceil_div(bd_square - total, denom)
        if trial > bd:
            bd, witness = trial, u
    return dict(L=length, n=n_taps, B_ph=pig, B_CS=cs, B_int=inte, B_bd=bd, u_bd=witness)

def delta(mask: int, length: int, stop_above: int | None = None) -> int:
    high = 0
    for d in range(1, length):
        high = max(high, (mask & (mask >> d)).bit_count())
        if stop_above is not None and high > stop_above:
            return high
    return high

def positions(mask: int) -> list[int]:
    return [i for i in range(mask.bit_length()) if (mask >> i) & 1]

def reflection(mask: int) -> int:
    span = mask.bit_length() - 1
    return sum(1 << (span - p) for p in positions(mask))

def exact_for_length(length: int) -> list[dict]:
    records = {n: dict(best=length, count=0, canonical=set(), witness=0) for n in range(2, length)}
    for right in range(1, 1 << (length - 1)):
        mask = 1 | (right << 1)
        n_taps = mask.bit_count()
        if n_taps == length:
            continue
        record = records[n_taps]
        value = delta(mask, length, record["best"])
        if value > record["best"]:
            continue
        if value < record["best"]:
            record.update(best=value, count=0, canonical=set(), witness=mask)
        record["count"] += 1
        record["canonical"].add(min(mask, reflection(mask)))
        if positions(mask) < positions(record["witness"]):
            record["witness"] = mask
    out = []
    for n_taps, item in records.items():
        b = bounds(length, n_taps)
        assert b["B_bd"] <= item["best"]
        out.append(b | dict(
            Delta=item["best"], gap=item["best"] - b["B_bd"],
            optimal_anchored_sets=item["count"],
            optimal_translation_reflection_classes=len(item["canonical"]),
            witness=positions(item["witness"]),
        ))
    return out
