"""Complement-capacity strengthening of the boundary relaxation.

The boundary relaxation of Theorem 1 constrains the window counts only from above by the
physical capacities h_x.  A tap set leaves exactly L-n positions empty, so every window
count is constrained from below as well:

        a_x = |T & W_x| >= h_x - (L-n)   for every window x.

Adding these lower bounds to the integer capacity problem can only increase the minimum
square sum, hence it yields a valid relaxation that is never weaker than B_bd:

        B_bd2(L,n) = max_u ceil( (Q_bd2(u) - n u) / (u(u-1)) )  >=  B_bd(L,n),

with Q_bd2 the minimum of sum_x a_x^2 over integers a with
        max(0, h_x - (L-n)) <= a_x <= h_x   and   sum_x a_x = n u.

Both Q_bd and Q_bd2 are computed by layer filling: the marginal cost of raising an entry
from level j to j+1 is 2j+1 independently of x, and the number of entries eligible for
that increment is  #{x : l_x <= j < h_x}.
"""

from __future__ import annotations

import json

import sys

def H(v, N, u):
    """#{x : h_x <= v}; the tent profile takes 1..u-1 twice and u exactly L-u+1 times."""
    if v <= 0:
        return 0
    if v >= u:
        return N
    return 2 * v

def profile_square_and_sum(f, N, u):
    """sum over windows of f(h_x)."""
    total = 0
    prev = 0
    for v in range(1, u + 1):
        cnt = H(v, N, u) - H(v - 1, N, u)
        total += cnt * f(v)
        prev = v
    return total

def Q_bd2(L, u, M, complement=None):
    c = (L - (M // u)) if False else None
    return None

def capacity_min(L, u, M, c):
    """Minimum square sum with lower bounds l_x = max(0, h_x - c) and upper bounds h_x."""
    N = L + u - 1
    # sum of l and sum of l^2
    sum_l = 0
    sum_l2 = 0
    if c < u:
        for v in range(c + 1, u + 1):
            cnt = H(v, N, u) - H(v - 1, N, u)
            sum_l += cnt * (v - c)
            sum_l2 += cnt * (v - c) * (v - c)
    R = M - sum_l
    if R < 0:
        return None
    total = sum_l2
    for j in range(0, u):
        avail = H(j + c, N, u) - H(j, N, u)
        if avail <= 0:
            continue
        take = min(avail, R)
        total += take * (2 * j + 1)
        R -= take
        if R == 0:
            break
    if R > 0:
        return None
    return total

def bound(L, n, use_complement=True):
    best = 0
    best_u = None
    for u in range(2, L + 1):
        M = n * u
        c = (L - n) if use_complement else u   # c >= u gives l_x = 0, i.e. Theorem 1
        Q = capacity_min(L, u, M, c)
        if Q is None:
            continue
        D = u * (u - 1)
        val = Q - M
        if val > 0:
            ceil = -((-val) // D)
            if ceil > best:
                best, best_u = ceil, u
    return best, best_u
