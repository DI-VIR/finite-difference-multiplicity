"""Declared prime-field modular seed variants with matched direct-delta refinement.

This is a witness experiment, not an exhaustive optimization of construction families.
Field polynomial arithmetic and irreducibility tests use SymPy.
"""
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
import json
import math
import random
import statistics
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from extremal import delta
from paley_construction import generate as paley_generate, subset
from sympy import divisors, factorint, isprime, primitive_root
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_irreducible_p, gf_pow_mod, gf_mul, gf_rem, gf_add

HERE = Path(__file__).resolve().parents[2]
PAIRS = [(24,11),(32,15),(40,19),(48,22),(64,29),(80,37),(160,71),(256,118),(320,142)]
FAMILIES = ('Paley', 'Singer', 'Ruzsa-Spence', 'RDS', 'Random')


@lru_cache(None)
def primitive_polynomial(q, degree):
    order = q**degree - 1
    factors = list(factorint(order))
    for coeffs in product(range(q), repeat=degree):
        if coeffs[-1] == 0:
            continue
        f = [1] + list(coeffs)
        if gf_irreducible_p(f, q, ZZ) and all(gf_pow_mod([1,0], order//p, f, q, ZZ) != [1] for p in factors):
            return f
    raise RuntimeError((q, degree))


@lru_cache(None)
def trace_circle(q, degree, target):
    f = primitive_polynomial(q, degree)
    basis_trace = []
    for j in range(degree):
        a, total = [1] + [0]*j, []
        for _ in range(degree):
            total = gf_add(total, a, q, ZZ)
            a = gf_pow_mod(a, q, f, q, ZZ)
        assert len(total) <= 1
        basis_trace.append(total[0] if total else 0)
    modulus = (q**3-1)//(q-1) if degree == 3 else q*q-1
    marks, a = [], [1]
    for i in range(modulus):
        tr = sum(c*basis_trace[len(a)-j-1] for j,c in enumerate(a)) % q
        if tr == target:
            marks.append(i)
        a = gf_rem(gf_mul(a, [1,0], q, ZZ), f, q, ZZ)
    assert len(marks) == (q+1 if degree == 3 else q)
    return tuple(marks), modulus, tuple(f)


def ruzsa_circle(p):
    theta = int(primitive_root(p))
    modulus = p*(p-1)
    # CRT: x=i mod p-1 and x=theta^i mod p.
    marks = sorted(i + (p-1)*((i-pow(theta,i,p)) % p) for i in range(p-1))
    return marks, modulus, theta


def modular_counts(marks, modulus):
    return Counter((b-a) % modulus for a in marks for b in marks if a != b)


def cut_and_thin(marks, modulus, n):
    marks = sorted(set(marks))
    gaps = [(marks[(i+1)%len(marks)]-a) % modulus for i,a in enumerate(marks)]
    largest = max(gaps)
    choices = [sorted((a-marks[(i+1)%len(marks)]) % modulus for a in marks)
               for i,gap in enumerate(gaps) if gap == largest]
    full = min(choices)
    return subset(full, n), largest


def select_seed(L, n, family):
    if family == 'Paley':
        p,g,gap,taps = paley_generate(L,n)
        original = sorted({x*x % p for x in range(1,p)})
        observed = max(modular_counts(original,p).values())
        assert observed <= g
        return taps, dict(q=p, original_modulus=p, modulus=p, fold_divisor=1,
            modular_cap=g, largest_gap=gap, original_marks=original,
            original_cardinality=len(original), projected_cardinality=len(original),
            discarded_collisions=0, guaranteed_cardinality=len(original),
            cardinality_check='theorem_and_explicit', actual_modular_max=observed)
    for q in range(max(3,n-1), 2*L+1):
        if not isprime(q):
            continue
        N = q*q+q+1 if family == 'Singer' else q*(q-1) if family == 'Ruzsa-Spence' else q*q-1
        count = q+1 if family == 'Singer' else q-1 if family == 'Ruzsa-Spence' else q
        for m in sorted((int(d) for d in divisors(N) if n <= d <= L), reverse=True):
            g = N//m
            guaranteed = count-g//2 if family != 'RDS' else count-(g-math.gcd(g,q-1))//2
            if guaranteed < n and family != 'Singer':
                continue
            if family == 'Singer':
                original,N,poly = trace_circle(q,3,0)
                details = dict(polynomial=list(poly), trace_target=0, field_degree=3)
            elif family == 'RDS':
                original,N,poly = trace_circle(q,2,1)
                details = dict(polynomial=list(poly), trace_target=1, field_degree=2)
            else:
                original,N,theta = ruzsa_circle(q)
                details = dict(primitive_root=theta)
            counts = modular_counts(original,N)
            assert max(counts.values()) == 1
            projected = sorted({a % m for a in original})
            assert len(projected) >= guaranteed
            if len(projected) < n:
                continue
            projected_counts = modular_counts(projected,m)
            observed = max(projected_counts.values())
            assert observed <= g
            taps,gap = cut_and_thin(projected,m,n)
            assert len(set(taps)) == n and min(taps) == 0 and max(taps) < L
            assert delta(sum(1<<t for t in taps),L) <= observed
            return taps, dict(q=q, original_modulus=N, modulus=m, fold_divisor=g,
                modular_cap=g, actual_modular_max=observed, original_cardinality=count,
                projected_cardinality=len(projected), discarded_collisions=count-len(projected),
                guaranteed_cardinality=guaranteed, cardinality_check='explicit' if guaranteed < n else 'theorem_and_explicit',
                largest_gap=gap, original_marks=list(original), **details)
    return None, dict(status='NOT_APPLICABLE_IN_DECLARED_PRIME_SCAN', q_limit=2*L)


def refine(L,n,initial,seed,budget):
    rng = random.Random(seed)
    current = list(initial)
    empty = sorted(set(range(L))-set(current))
    bits = sum(1<<t for t in current)
    value = delta(bits,L)
    for _ in range(budget):
        i,j = rng.randrange(n),rng.randrange(L-n)
        x,y = current[i],empty[j]
        candidate = bits^(1<<x)^(1<<y)
        score = delta(candidate,L,stop_above=value)
        if score <= value:
            bits,value = candidate,score
            current[i],empty[j] = y,x
    assert value == delta(bits,L)
    return value, sorted(current)


def run(seeds=30,budget=3000):
    started = time.perf_counter()
    records, summary, constructors = [], [], []
    prior = json.loads((HERE/'data/literature_comparison.json').read_text())
    lower = {(r['L'],r['n']):max(r['B_bd'],r['table_lower']) for r in prior}
    for L,n in PAIRS:
        for family in FAMILIES:
            initial,details = (None,{}) if family == 'Random' else select_seed(L,n,family)
            constructors.append(dict(L=L,n=n,family=family,initial_taps=initial,**details))
            if initial is None and family != 'Random':
                print('NOT APPLICABLE',L,n,family,flush=True)
                continue
            for s in range(seeds):
                seed = 2026100200+s
                taps = random.Random(seed).sample(range(L),n) if family == 'Random' else initial
                before = delta(sum(1<<t for t in taps),L)
                after,witness = refine(L,n,taps,seed,budget)
                independent = Counter(b-a for i,a in enumerate(witness) for b in witness[i+1:])
                assert max(independent.values()) == after and after >= lower[(L,n)]
                records.append(dict(L=L,n=n,family=family,seed=seed,budget=budget,
                    before=before,after=after,initial_taps=sorted(taps),taps=witness,B_comb=lower[(L,n)]))
            cell = [r for r in records if (r['L'],r['n'],r['family']) == (L,n,family)]
            values = [r['after'] for r in cell]
            summary.append(dict(L=L,n=n,family=family,B_comb=lower[(L,n)],seeds=seeds,budget=budget,
                before_mean=statistics.mean(r['before'] for r in cell),after_mean=statistics.mean(values),
                after_SD=statistics.stdev(values),best=min(values),remaining_gap=min(values)-lower[(L,n)]))
            print('DONE',L,n,family,'best',min(values),'mean',round(statistics.mean(values),3),flush=True)
            for name,data in [('modular_runs',records),('modular_summary',summary),('modular_constructors',constructors)]:
                (HERE/'data'/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    report = dict(status='PASS',runs=len(records),cells=len(summary),seeds=seeds,budget=budget,
        seconds=round(time.perf_counter()-started,3),source='https://arxiv.org/html/2607.07931v1',
        scope='Declared prime-field variants only; not an exhaustive or optimal parameter search over all construction families.')
    (HERE/'tests/modular_comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seeds',type=int,default=30)
    parser.add_argument('--budget',type=int,default=3000)
    args = parser.parse_args()
    run(args.seeds,args.budget)
