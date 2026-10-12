"""Exact final-closure checks; no optimization result or human review is fabricated."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path
import argparse
import json
import time

from extremal import bounds, ceil_div, min_boundary_square
from constructions.modular_comparison import (
    cut_and_thin, modular_counts, select_seed, trace_circle, ruzsa_circle,
)
from sympy import factorint, isprime
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_irreducible_p, gf_pow_mod

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    return json.loads((ROOT / 'data' / (name + '.json')).read_text(encoding='utf-8'))


def audit(complete_metadata=False):
    started = time.perf_counter()
    constructors = load('modular_constructors')
    if complete_metadata:
        for row in constructors:
            if row['family'] == 'Paley':
                taps, details = select_seed(row['L'], row['n'], row['family'])
                assert taps == row['initial_taps']
                assert all(details[k] == v for k, v in row.items() if k in details)
                row.update(details)
        (ROOT / 'data/modular_constructors.json').write_text(
            json.dumps(constructors, indent=2) + '\n', encoding='utf-8')

    boundary_windows = 0
    for L in range(2,25):
        for n in range(2,L+1):
            correction_exceeds_slack = False
            integer_bound = bounds(L,n)['B_int']
            for u in range(2,L+1):
                N, M = L+u-1, n*u
                caps = [max(0,min(L,x+u)-max(0,x)) for x in range(1-u,L)]
                assert caps == list(range(1,u))+[u]*(L-u+1)+list(range(u-1,0,-1))
                assert len(caps) == N and sum(caps) == L*u
                Q, k, eta = min_boundary_square(L,n,u)
                assert k*(L+u-k) <= M
                if k < u:
                    assert 0 <= eta < L+u-2*k-1
                else:
                    assert eta == 0 and M == L*u
                profile = [min(h,k) for h in caps]
                eligible = [i for i,h in enumerate(caps) if h >= k+1]
                assert eta <= len(eligible)
                for i in eligible[:eta]:
                    profile[i] += 1
                assert sum(profile) == M and sum(x*x for x in profile) == Q
                q, r = divmod(M,N)
                Qint = N*q*q+(2*q+1)*r
                assert (Q>Qint) == (q>=2 or q==1 and r>N-2)
                sigma = integer_bound*u*(u-1)-(Qint-M)
                assert sigma >= 0
                correction_exceeds_slack |= Q-Qint > sigma
                boundary_windows += 1
            assert (bounds(L,n)['B_bd']>integer_bound) == correction_exceeds_slack

    examples = []
    for L,n,u,qi,qb,ci,cb in [(8,2,3,6,6,0,0),(7,4,3,18,18,1,1),
                              (6,5,3,29,31,3,3),(8,7,4,74,80,4,5)]:
        N,M = L+u-1,n*u
        q,r = divmod(M,N)
        assert N*q*q+(2*q+1)*r == qi
        assert min_boundary_square(L,n,u)[0] == qb
        assert ceil_div(qi-M,u*(u-1)) == ci and ceil_div(qb-M,u*(u-1)) == cb
        examples.append(dict(L=L,n=n,u=u,q=q,eta=r,Q_int=qi,Q_bd=qb,
                             local_integer_ceiling=ci,local_boundary_ceiling=cb))

    slack_sets = 0
    expectation_totals = Counter()
    for L in range(2,9):
        for n in range(2,L+1):
            for T in combinations(range(L),n):
                multiplicity = Counter(b-a for i,a in enumerate(T) for b in T[i+1:])
                actual_cap = max(multiplicity.values())
                for u in range(2,L+1):
                    a = [sum(x<=t<x+u for t in T) for x in range(1-u,L)]
                    expectation_totals[L,n,u] += sum(x*x for x in a)
                    Q = min_boundary_square(L,n,u)[0]
                    D,M = u*(u-1),n*u
                    g = ceil_div(Q-M,D)
                    assert sum(a) == M
                    assert sum(x*x for x in a)-M == 2*sum((u-d)*multiplicity[d] for d in range(1,u))
                    if g == actual_cap:
                        s = D*g-(Q-M)
                        loss = sum(x*x for x in a)-Q
                        deficit = 2*sum((u-d)*(g-multiplicity[d]) for d in range(1,u))
                        assert 0 <= s < D and loss>=0 and deficit>=0 and loss+deficit==s
                        if s == 0:
                            assert loss==0 and all(multiplicity[d]==g for d in range(1,u))
                slack_sets += 1
    for (L,n,u),total in expectation_totals.items():
        p=Fraction(n*(n-1),L*(L-1))
        expected=n*u+p*u*(u-1)*(L-Fraction(u+1,3))
        assert Fraction(total,comb(L,n))==expected

    family = load('infinite_family_checks')
    assert [r['k'] for r in family] == list(range(3,41))+[50,64,80,100]
    for row in family:
        k,L,n = row['k'],row['L'],row['n']
        assert (L,n)==(4*k*k,2*k*k) and Fraction(n,L)==Fraction(1,2)
        b = bounds(L,n)
        assert (b['B_int'],b['B_bd'])==(row['B_int'],row['B_bd'])
        Q,layer,eta = min_boundary_square(L,n,2*k)
        assert (layer,eta)==(k-1,3*k*k+1)
        assert Q==Fraction(12*k**4-4*k**3+3*k*k+k,3)
        assert Fraction(Q-n*2*k,2*k*(2*k-1))==k*k-Fraction(5*k+1,6)
        u=2*k+1;N=L+u-1;M=n*u
        assert Fraction(M*M-M*N,N*u*(u-1))==k*k-k
        assert ceil_div(k-1,6)-2 <= row['actual_gain'] <= k
        if k>=14:
            assert row['actual_gain']>0

    construction_rows = 0
    eligibility = Counter()
    for row in constructors:
        if row['family']=='Random':
            assert row['initial_taps'] is None
            continue
        q,N,m = row['q'],row['original_modulus'],row['modulus']
        original = row['original_marks']
        assert isprime(q) and N==m*row['fold_divisor']
        if row['family'] in ('Singer','RDS'):
            f,degree = row['polynomial'],row['field_degree']
            assert len(f)==degree+1 and gf_irreducible_p(f,q,ZZ)
            order=q**degree-1
            assert gf_pow_mod([1,0],order,f,q,ZZ)==[1]
            assert all(gf_pow_mod([1,0],order//p,f,q,ZZ)!=[1] for p in factorint(order))
            regenerated,modulus,polynomial = trace_circle(q,degree,row['trace_target'])
            assert (list(regenerated),modulus,list(polynomial))==(original,N,f)
        elif row['family']=='Ruzsa-Spence':
            regenerated,modulus,theta = ruzsa_circle(q)
            assert (regenerated,modulus,theta)==(original,N,row['primitive_root'])
        else:
            assert original == sorted({x*x % q for x in range(1,q)})
        original_max = max(modular_counts(original,N).values())
        if row['family']=='Paley':
            assert original_max <= row['modular_cap']
        else:
            assert original_max == 1
        projected = sorted({a % m for a in original})
        assert len(original)==row['original_cardinality']
        assert len(projected)==row['projected_cardinality']>=row['n']
        assert len(original)-len(projected)==row['discarded_collisions']
        assert len(projected)>=row['guaranteed_cardinality']
        if row['cardinality_check']=='explicit':
            assert row['family']=='Singer' and row['guaranteed_cardinality']<row['n']
        else:
            assert row['guaranteed_cardinality']>=row['n']
        eligibility[row['cardinality_check']] += 1
        observed = max(modular_counts(projected,m).values())
        assert observed == row['actual_modular_max'] <= row['modular_cap']
        taps,gap = cut_and_thin(projected,m,row['n'])
        assert (taps,gap)==(row['initial_taps'],row['largest_gap'])
        assert len(set(taps))==row['n'] and min(taps)==0 and max(taps)<row['L']
        construction_rows += 1
    assert construction_rows==36
    runs = load('modular_runs')
    summary = load('modular_summary')
    assert len(runs)==1350 and len(summary)==45
    for cell in summary:
        values = [r for r in runs if (r['L'],r['n'],r['family'])==(cell['L'],cell['n'],cell['family'])]
        assert len(values)==30 and {r['seed'] for r in values}==set(range(2026100200,2026100230))
        assert {r['budget'] for r in values}=={3000}
        assert min(r['after'] for r in values)==cell['best']

    from regenerate_tables import read_tables
    catalog = read_tables()
    if not complete_metadata:
        assert catalog['4']['headers'][7]=='Closed by B_real?'
        assert len(catalog['S14']['rows'])==23
        assert [int(r[0]) for r in catalog['S14']['rows'][:20]]==list(range(3,23))
        milp = load('milp_complete')
        assert len(catalog['6']['rows'])==len(milp)==4
        for record,row in zip(milp,catalog['6']['rows']):
            g = bounds(record['L'],record['n'])['B_bd']
            assert record['target']==g==record['exact_Delta'] and row[2]==f'B_bd = {g}'
    grid = load('bound_landscape')
    assert all(r['n']==round(.45*r['L']) for r in grid if r['target_density']==.45)
    for L in range(4,65):
        T=[x for x in range(L) if x!=1]
        assert max(Counter(b-a for i,a in enumerate(T) for b in T[i+1:]).values())==L-3
        if L>=9:
            T=[x for x in range(L) if x not in (2,L-3)]
            assert max(Counter(b-a for i,a in enumerate(T) for b in T[i+1:]).values())==L-5
    return dict(status='PASS',boundary_windows=boundary_windows,slack_identity_sets=slack_sets,
        expectation_domains=len(expectation_totals),
        hand_examples=examples,infinite_family_checks=len(family),algebraic_constructor_records=construction_rows,
        construction_eligibility=dict(eligibility),construction_cells=len(summary),construction_runs=len(runs),
        medium_solver_records=4,metadata_only_completion=complete_metadata,
        human_proof_review='NOT_PERFORMED_BY_THIS_AUTOMATED_CHECK',
        seconds=round(time.perf_counter()-started,3))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--complete-constructor-metadata',action='store_true')
    args=parser.parse_args()
    report=audit(args.complete_constructor_metadata)
    (ROOT/'tests/final_details_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
