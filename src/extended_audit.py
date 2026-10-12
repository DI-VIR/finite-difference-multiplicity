"""Independent checks of the new revision records without rerunning searches."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import statistics
import time
from binary_window_milp import square_sum
from extremal import bounds, min_boundary_square
from fourier_verify import run as fourier_verify

ROOT=Path(__file__).resolve().parents[1]


def load(name): return json.loads((ROOT/'data'/(name+'.json')).read_text())


def audit():
    start=time.perf_counter()
    fourier=fourier_verify()
    runs=load('seed_budget_runs')
    assert len(runs)==4860
    keys={(r['L'],r['n'],r['family'],r['seed'],r['budget']):r for r in runs}
    assert len(keys)==len(runs)
    for r in runs:
        T=r['witness'];L,n=r['L'],r['n']
        assert len(T)==len(set(T))==n and T==sorted(T) and 0<=min(T)<=max(T)<L
        assert max(Counter(b-a for i,a in enumerate(T) for b in T[i+1:]).values())==r['value']>=bounds(L,n)['B_bd']
        if r['budget']==0: assert T==r['initial']
    old=load('modular_runs')
    for r in old:
        new=keys[(r['L'],r['n'],r['family'],r['seed'],3000)]
        assert new['value']==r['after'] and new['witness']==r['taps']
    for r in load('seed_budget_summary'):
        cell=[v['value'] for v in runs if (v['L'],v['n'],v['family'],v['budget'])==(r['L'],r['n'],r['family'],r['budget'])]
        assert len(cell)==30 and r['mean']==statistics.mean(cell) and r['SD']==statistics.stdev(cell) and r['best']==min(cell)
    for r in load('seed_budget_contrasts'):
        d=[keys[(r['L'],r['n'],r['family'],seed,r['budget'])]['value']-
           keys[(r['L'],r['n'],'Random',seed,r['budget'])]['value'] for seed in range(2026100200,2026100230)]
        assert r['algebraic_minus_random_mean']==statistics.mean(d)
        assert r['paired_bootstrap_95_low']<=r['paired_bootstrap_95_high']
    for r in load('window_scalability'):
        assert len(r['witness'])==len(set(r['witness']))==r['n']
        q=square_sum(r['witness'],r['u']);qbd=min_boundary_square(r['L'],r['n'],r['u'])[0]
        assert q==r['Q_upper']>=qbd and r['time_limit']==60 and r['L']>24
        if r['status']==0: assert q==qbd
        else: assert r['status']==1 and r['outcome']=='TIME_LIMIT_BRACKET'
    medium={(r['L'],r['n']):r for r in load('verified_medium')}
    for r in load('medium_witness_provenance'):
        archived=medium[(r['L'],r['n'])]
        assert r['witness']==list(map(int,archived['witness'].split())) and r['value']==archived['best_feasible']
        assert r['witness_sha256']==sha256(json.dumps(r['witness'],separators=(',',':')).encode()).hexdigest()
        assert r['source_record_index']>=0 and len(r['original_file_sha256'])==64
    for r in load('xu_comparison'):
        L,n=r['L'],r['n']
        values=[Fraction(0)]
        for D in range(2,L+1):
            values.extend([Fraction(n*n,L+D)-Fraction(n+1,D),Fraction(2*(D*n-L),D*(D-1))-Fraction(2,D)])
        assert r['B_Xu']==max(-(-v.numerator//v.denominator) for v in values)
    report=dict(status='PASS',seconds=round(time.perf_counter()-start,3),
        exact_Fourier_cases=fourier['cases'],exact_Fourier_leaves=sum(r['exact_leaves'] for r in fourier['records']),
        new_checkpoint_witnesses=len(runs),unchanged_3000_prefix_witnesses=len(old),
        larger_binary_window_runs=6,analytically_closed_window_minima=2,time_limit_brackets=4,
        historical_coordinate_source_matches=9,
        scope='Saved certificates and populations independently checked; human proof review still pending.')
    (ROOT/'tests/extended_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return report


if __name__=='__main__': audit()
