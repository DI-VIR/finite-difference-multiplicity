"""Audit of the finite-table comparison: corrected counts, reverse cases, combined certificate.

The manuscript previously stated 379 / 83 / 0 for
(B_bd > prior) / (equal) / (reverse).  The archived CSV gives 379 / 74 / 9.  This module
recomputes the counts from the raw grid, stratifies them by coverage regime, audits every
reverse case back to the source table cell, and builds the combined certificate

    B_comb(L,n) = max(B_bd(L,n), B_prior(L,n)).

Exhaustive anchored enumeration is used for the small grid rows where it is affordable,
so that B_bd, B_prior, B_comb can be compared against the true value of Delta(L,n).
"""

from __future__ import annotations

import itertools

import json

import sys

from collections import Counter

from common import HERE, save, csv_read

def delta_of(taps, L):
    bits = 0
    for t in taps:
        bits |= 1 << t
    best = 0
    for d in range(1, L):
        v = (bits & (bits >> d)).bit_count()
        if v > best:
            best = v
            if best == len(taps) - 1:
                break
    return best

def exact_delta(L, n, lower, cap=6_000_000):
    """Enumerate anchored n-subsets; stop as soon as the lower bound is attained."""
    best = None
    count = 0
    for combo in itertools.combinations(range(1, L), n - 1):
        taps = (0,) + combo
        d = delta_of(taps, L)
        count += 1
        if best is None or d < best:
            best = d
            if best <= lower:
                return best, count, True
        if count >= cap:
            return best, count, False
    return best, count, True

def source_entries(L, n, table):
    """For each excluded g, the tabulated cell (g,k,diameter) with smallest k that excludes it."""
    out = []
    for g in sorted({row["g"] for row in table}):
        hits = [(row["k"], row["diameter"]) for row in table
                if row["g"] == g and row["k"] <= n and row["diameter"] > L - 1]
        if hits:
            k, diam = min(hits)
            out.append(dict(g=g, k=k, diameter=diam, margin=diam - (L - 1)))
    return out
