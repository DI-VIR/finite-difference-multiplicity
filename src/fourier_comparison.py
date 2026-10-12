"""Fourier cuts for the declared integer multiplicity relaxation.

Negative polynomial evaluations are exact rational certificates. MILP
infeasibility is a numerical solver outcome, not an exact Farkas certificate.
Finite Toeplitz PSD, Fourier nonnegativity and Boolean moment SDP are distinct.
"""
from fractions import Fraction
from pathlib import Path
import json
import time

import numpy as np
from numpy.polynomial import chebyshev as ch
from scipy.linalg import toeplitz
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from sympy import Poly, Rational, Symbol, chebyshevt

ROOT = Path(__file__).resolve().parents[1]
X = Symbol('x')


def coefficients(x, degree):
    x = Fraction(x)
    seq = [Fraction(1), x]
    for _ in range(2, degree + 1):
        seq.append(2*x*seq[-1]-seq[-2])
    return [2*t for t in seq[1:degree+1]]


def exact_value(n, r, x):
    return Fraction(n) + sum(a*b for a, b in zip(r, coefficients(x, len(r))))


def separate(n, r):
    c = np.r_[float(n), 2*np.asarray(r, dtype=float)]
    roots = ch.chebroots(ch.chebder(c))
    xs = [-1.0, 1.0] + [float(z.real) for z in roots
                              if abs(z.imag) < 1e-7 and -1 < z.real < 1]
    x = min(xs, key=lambda v: ch.chebval(v, c))
    rational = Fraction(x).limit_denominator(1000000)
    if exact_value(n, r, rational) < 0:
        return rational, False
    # Root isolation of the integer polynomial checks signs between all roots.
    p = Poly(n + 2*sum(a*chebyshevt(i, X) for i, a in enumerate(r, 1)), X)
    intervals = p.intervals(eps=Rational(1, 10**8))
    points, last = [Rational(-1)], Rational(-1)
    for (a, b), _ in intervals:
        if b < -1 or a > 1:
            continue
        a, b = max(a, Rational(-1)), min(b, Rational(1))
        if a > last:
            points.append((last+a)/2)
        last = b
    if last < 1:
        points.append((last+1)/2)
    points.append(Rational(1))
    for point in points:
        if p.eval(point) < 0:
            return Fraction(int(point.p), int(point.q)), False
    return None, True


def rational_exclusion(base, lower, upper, cuts, n, limits, mass):
    """Validate a nonnegative rational inequality combination with bounded residuals."""
    degree = len(limits)
    A, b = [], []
    for row, lo, hi in zip(base, lower, upper):
        row = [Fraction(int(v)) for v in row]
        if np.isfinite(lo):
            A.append([-v for v in row]); b.append(Fraction(-int(lo)))
        if np.isfinite(hi):
            A.append(row); b.append(Fraction(int(hi)))
    odd = [Fraction(int(d%2==1)) for d in range(1,degree+1)]
    A.extend([odd,[-v for v in odd]]); b.extend([Fraction(mass),Fraction(-mass)])
    for x in cuts:
        A.append([-v for v in coefficients(x,degree)]); b.append(Fraction(n))
    for i, limit in enumerate(limits):
        row = [Fraction(int(j==i)) for j in range(degree)]
        A.append(row); b.append(Fraction(limit))
        A.append([-v for v in row]); b.append(Fraction(0))
    matrix = np.asarray([[float(v) for v in row] for row in A])
    rhs = np.asarray([float(v) for v in b])
    result = linprog(rhs, A_eq=np.vstack([matrix.T,np.ones(len(b))]),
        b_eq=np.r_[np.zeros(degree),1.0], bounds=(0,None),method='highs')
    if not result.success or result.fun >= -1e-8:
        return dict(outcome='NO_LP_EXCLUSION',message=result.message)
    weights = [Fraction(float(max(0,v))).limit_denominator(10**9) for v in result.x]
    active = [(i,w) for i,w in enumerate(weights) if w]
    combined = [sum(w*A[i][j] for i,w in active) for j in range(degree)]
    total = sum(w*b[i] for i,w in active)
    floor = sum(min(v*limit,Fraction(0)) for v,limit in zip(combined,limits))
    if total >= floor:
        return dict(outcome='RATIONAL_VALIDATION_FAILED')
    return dict(outcome='EXACT_RATIONAL_EXCLUSION',
        weights=[[i,str(w)] for i,w in active],
        combined_rhs=str(total), residual_lower_bound=str(floor),
        strictly_positive_margin=str(floor-total),
        cuts=[str(x) for x in cuts],
        proof='A r <= b implies sum(y_i b_i) >= min over 0<=r<=limits of (y A) r; the opposite strict inequality is checked in Fraction arithmetic.')


def solve_case(case, profiles, time_limit=60):
    L, n, g = case['L'], case['n'], case['cap']
    degree = L-1
    base = [np.ones(degree)]
    lower = [n*(n-1)/2]
    upper = lower.copy()
    for row in profiles:
        u = row['u']
        base.append(np.asarray([2*max(0, u-d) for d in range(1, L)], dtype=float))
        lower.append(row['Q_real']-n*u)
        upper.append(np.inf)
    odd = np.asarray([int(d % 2 == 1) for d in range(1, L)], dtype=float)
    masses = sorted({e*(n-e) for e in range(max(0, n-L//2), min(n, (L+1)//2)+1)})
    outcomes = []
    start = time.perf_counter()
    for mass in masses:
        cuts = [Fraction(i, 128) for i in range(-128, 129)]
        attempts = []
        for iteration in range(100):
            matrix = [*base, odd, *[[float(v) for v in coefficients(x, degree)] for x in cuts]]
            lo = [*lower, mass, *([-n]*len(cuts))]
            hi = [*upper, mass, *([np.inf]*len(cuts))]
            tick = time.perf_counter()
            result = milp(np.zeros(degree), integrality=np.ones(degree),
                bounds=Bounds(np.zeros(degree), [min(g, n, L-d) for d in range(1,L)]),
                constraints=LinearConstraint(np.asarray(matrix), lo, hi),
                options={'time_limit':time_limit, 'mip_rel_gap':0.0})
            attempt = dict(iteration=iteration, status=int(result.status),
                message=result.message, seconds=round(time.perf_counter()-tick,6), cuts=len(cuts))
            attempts.append(attempt)
            if result.status != 0:
                certificate = rational_exclusion(base,lower,upper,cuts,n,
                    [min(g,n,L-d) for d in range(1,L)],mass) if result.status==2 else None
                outcomes.append(dict(odd_mass=mass, attempts=attempts,
                    outcome='EXACT_RATIONAL_EXCLUSION' if certificate and certificate['outcome']=='EXACT_RATIONAL_EXCLUSION' else 'NUMERICAL_INFEASIBLE' if result.status == 2 else 'UNRESOLVED',
                    rational_certificate=certificate))
                break
            r = np.rint(result.x).astype(int).tolist()
            assert np.allclose(result.x, r, atol=1e-6, rtol=0)
            assert sum(r)==n*(n-1)//2 and sum(r[::2])==mass
            assert all(0 <= v <= min(g,n,L-d) for d,v in enumerate(r,1))
            assert all(sum(2*max(0,p['u']-d)*v for d,v in enumerate(r,1)) >= p['Q_real']-n*p['u'] for p in profiles)
            x, certified = separate(n,r)
            if certified:
                outcomes.append(dict(odd_mass=mass, attempts=attempts,
                    outcome='EXACT_FOURIER_FEASIBLE', witness=r,
                    proof='Integer Chebyshev polynomial real-root isolation and rational signs'))
                break
            attempt['negative_cut'] = str(x)
            attempt['exact_negative_value'] = str(exact_value(n,r,x))
            assert x not in cuts
            cuts.append(x)
        else:
            outcomes.append(dict(odd_mass=mass, attempts=attempts, outcome='ITERATION_LIMIT'))
    feasible = any(r['outcome']=='EXACT_FOURIER_FEASIBLE' for r in outcomes)
    all_exact = all(r['outcome']=='EXACT_RATIONAL_EXCLUSION' for r in outcomes)
    all_infeasible = all(r['outcome'] in ('NUMERICAL_INFEASIBLE','EXACT_RATIONAL_EXCLUSION') for r in outcomes)
    return dict(L=L,n=n,cap=g,parity_branches=outcomes,
        outcome='EXACT_FOURIER_FEASIBLE' if feasible else 'EXACT_RATIONAL_EXCLUSION' if all_exact else 'NUMERICAL_INFEASIBLE' if all_infeasible else 'UNRESOLVED',
        seconds=round(time.perf_counter()-start,6),
        interpretation='No exact integer infeasibility certificate is claimed.')


def run():
    cases = json.loads((ROOT/'data/residual_diagnostics.json').read_text())
    profiles = {(r['L'],r['n']):r['per_u'] for r in json.loads((ROOT/'data/seven_case_profiles.json').read_text())}
    test_xs = ['241/250','41/50','102/125','-241/250','443/500']
    records = []
    for case, x in zip(cases,test_xs):
        n,L = case['n'],case['L']
        r = case['strengthened_witness']
        value = exact_value(n,r,Fraction(x))
        assert value < 0
        record = dict(L=L,n=n,cap=case['cap'],cos_theta=x,
            exact_negative_value=str(value),negative_value=float(value),
            finite_Toeplitz_min_eigenvalue=float(np.linalg.eigvalsh(toeplitz([n,*r]))[0]),
            finite_Toeplitz_status='Numerical diagnostic only',
            new_relaxation=solve_case(case,profiles[(L,n)]))
        records.append(record)
        (ROOT/'data/fourier_comparison.json').write_text(json.dumps(records,indent=2)+'\n')
        print(L,n,record['new_relaxation']['outcome'],record['new_relaxation']['seconds'],flush=True)
    print('Complete; all five original-vector exclusions have exact rational sign certificates.',flush=True)


if __name__ == '__main__':
    run()
