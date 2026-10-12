"""Regression checks of new finite witnesses and the analytic family formulas."""
from pathlib import Path
from collections import Counter
from itertools import product
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from extremal import bounds, min_boundary_square
from realizability_bound import tent_capacity

# Exhaustively check the strict square-sum criterion in small capacitated domains.
for L in range(2,7):
    for u in range(2,L+1):
        h = tent_capacity(L,u).tolist()
        minimum = {}
        for a in product(*(range(x+1) for x in h)):
            mass,cost = sum(a),sum(x*x for x in a)
            minimum[mass] = min(cost,minimum.get(mass,cost))
        for n in range(2,L+1):
            N,M = L+u-1,n*u
            q,eta = divmod(M,N)
            qint = N*q*q+(2*q+1)*eta
            qbd = min_boundary_square(L,n,u)[0]
            assert qbd == minimum[M]
            assert (qbd>qint) == (q>=2 or q==1 and eta>N-2)

for r in json.loads((ROOT/'data/infinite_family_checks.json').read_text()):
    k = r['k'];b = bounds(4*k*k,2*k*k)
    assert r['actual_gain'] == b['B_bd']-b['B_int']
    assert r['proved_gain'] <= r['actual_gain'] <= k
for r in json.loads((ROOT/'data/residual_diagnostics.json').read_text()):
    L,n,g,w = r['L'],r['n'],r['cap'],r['strengthened_witness']
    assert len(w) == L-1 and sum(w) == n*(n-1)//2 and max(w) == g
    assert all(0<=x<=min(g,n,L-d) for d,x in enumerate(w,1))
    assert sum(w[::2]) in {e*(n-e) for e in range(n+1)}
for r in json.loads((ROOT/'data/modular_runs.json').read_text()):
    T = r['taps'];hist = Counter(b-a for i,a in enumerate(T) for b in T[i+1:])
    assert len(T)==len(set(T))==r['n'] and min(T)>=0 and max(T)<r['L']
    assert max(hist.values())==r['after']>=r['B_comb']
print('PASS: layer strictness, infinite-family instances, strengthened vectors and 1350 construction witnesses.')
