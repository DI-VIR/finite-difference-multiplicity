"""Independent integer checks, optional complete enumeration, and table export."""
from pathlib import Path
from collections import Counter
import argparse
import csv
import json
import time
from common import HERE
from extremal import bounds, delta
from realizability_bound import enumerate_L, B_from_Q, Q_bd, Q_walk
from np_reduction import verify as verify_reduction
from regenerate_tables import regenerate


def load(name):return json.loads((HERE/'data'/(name+'.json')).read_text())


def check(full=False):
    started=time.perf_counter();records=load('exact_map_nontrivial')
    binary_windows={(r['L'],r['n'],r['u']):r['Q_real'] for r in load('exact_binary_window_minima')}
    assert len(records)==252 and {(r['L'],r['n']) for r in records}=={(L,n) for L in range(4,25) for n in range(2,L)}
    for r in records:
        L,n=r['L'],r['n'];T=r['witness'];cc=Counter(b-a for i,a in enumerate(T) for b in T[i+1:])
        assert len(T)==len(set(T))==n and min(T)==0 and max(T)<L
        assert max(cc.values())==r['Delta']==delta(sum(1<<t for t in T),L)
        assert bounds(L,n)['B_bd']==r['B_bd']<=r['B_real']<=r['Delta']
    assert sum(r['B_bd']==r['Delta'] for r in records)==245
    assert sum(r['B_real']==r['Delta'] for r in records)==247
    trivial=load('exact_map_trivial')
    assert len(trivial)==21
    for r in trivial:
        L=r['L'];assert r['n']==L and r['Delta']==L-1 and r['witness']==list(range(L))
        assert bounds(L,L)['B_bd']==L-1
    audit=load('realizability_audit');residual=[r for r in audit if r.get('hierarchy_exhausted')]
    assert len(residual)==5
    for r in residual:
        vec=r['integer_r_at_B_bd'];assert len(vec)==r['L']-1 and min(vec)>=0
        assert max(vec)==r['B_bd']<r['Delta'] and sum(vec)==r['n']*(r['n']-1)//2
        assert all(2*sum((p['u']-d)*vec[d-1] for d in range(1,p['u']))>=p['Q_real']-r['n']*p['u'] for p in r['per_u'])
    for r in load('residual_diagnostics'):
        L,n,g,w=r['L'],r['n'],r['cap'],r['strengthened_witness']
        assert sum(w)==n*(n-1)//2 and max(w)==g and len(w)==L-1
        assert all(0<=a<=min(g,n,L-d) for d,a in enumerate(w,1))
        assert sum(w[::2]) in {e*(n-e) for e in range(n+1)}
        original=next(a for a in residual if (a['L'],a['n'])==(L,n))
        assert all(2*sum((p['u']-d)*w[d-1] for d in range(1,p['u']))>=p['Q_real']-n*p['u'] for p in original['per_u'])
    walk_windows=0
    for r in load('seven_case_profiles'):
        for p in r['per_u']:
            L,n,u=r['L'],r['n'],p['u']
            assert p['Q_bd']==Q_bd(L,u,n)<=p['Q_walk']==Q_walk(L,u,n)<=p['Q_real']==binary_windows[(L,n,u)]
            assert sum(p['walk_profile'])==n*u and sum(a*a for a in p['walk_profile'])==p['Q_walk']
            walk_windows+=1
    joint=load('joint_bound_audit');assert len(joint)==714
    for r in joint:
        L,n=r['L'],r['n'];v=r['integer_witness'];m=r['max_bd_ph']
        assert r['joint_status']=='feasible' and sum(v)==n*(n-1)//2 and max(v)<=m and min(v)>=0
        assert all(2*sum((u-d)*v[d-1] for d in range(1,u))>=Q_bd(L,u,n)-n*u for u in range(2,L+1))
    for r in load('milp_complete'):
        if r.get('witness'):
            T=r['witness'];assert len(set(T))==r['n'] and min(T)>=0 and max(T)<r['L']
            assert delta(sum(1<<t for t in T),r['L'])<=r['target']
        if r.get('exact_Delta') is not None:
            assert bounds(r['L'],r['n'])['B_bd']==r['exact_Delta']
    for r in load('paley_complete'):
        assert delta(sum(1<<t for t in r['initial_taps']),r['L'])==r['before']<=r['modular_g']
        assert delta(sum(1<<t for t in r['taps']),r['L'])==r['after']
    modular=load('modular_runs');assert len(modular)==1350
    assert len({(r['L'],r['n'],r['family']) for r in modular})==45
    for r in modular:
        T=r['taps'];assert len(T)==len(set(T))==r['n'] and min(T)>=0 and max(T)<r['L']
        assert max(Counter(b-a for i,a in enumerate(T) for b in T[i+1:]).values())==r['after']>=r['B_comb']
        assert delta(sum(1<<t for t in r['initial_taps']),r['L'])==r['before']
    grid=load('bound_landscape')
    assert len(grid)==462
    for r in grid:
        b=bounds(r['L'],r['n'])
        assert all(b[k]==r[k] for k in ('B_ph','B_CS','B_int','B_bd','u_bd'))
        assert r['bd_minus_int']==b['B_bd']-b['B_int']
    assert sum(r['bd_minus_int']>0 for r in grid)==189
    published=load('published_diameters');literature=load('literature_grid_comparison')
    splits={'direct':[0,0,0],'transferred':[0,0,0]}
    for r in literature:
        eligible=sorted({p['g'] for p in published if p['n']<=r['n'] and p['diameter']>r['L']-1})
        assert eligible==r['excluded_g']
        assert r['table_lower']==(max(eligible)+1 if eligible else 0)
        direct=any(p['n']==r['n'] for p in published)
        sign=0 if r['B_bd']>r['table_lower'] else 1 if r['B_bd']==r['table_lower'] else 2
        splits['direct' if direct else 'transferred'][sign]+=1
    assert splits=={'direct':[191,74,9],'transferred':[188,0,0]},splits
    for r in load('infinite_family_checks'):
        k=r['k'];b=bounds(r['L'],r['n'])
        assert b['B_bd']-b['B_int']==r['actual_gain'] and r['proved_gain']<=r['actual_gain']<=k
    with (HERE/'data/verified_medium.csv').open(encoding='utf-8-sig',newline='') as f:
        medium=list(csv.DictReader(f))
    assert len(medium)==9
    for r in medium:
        T=list(map(int,r['witness'].split()));L,n=int(r['L']),int(r['n'])
        assert len(T)==len(set(T))==n and min(T)>=0 and max(T)<L
        assert delta(sum(1<<t for t in T),L)==int(r['best_feasible'])
        assert bounds(L,n)['B_bd']==int(r['B_bd'])<=int(r['best_feasible'])
    brute_pairs=0;optimum_count_pairs=0
    for L in range(4,25 if full else 13):
        D,Q,evidence=enumerate_L(L,with_evidence=True)
        for r in (x for x in records if x['L']==L):
            n=r['n'];e=evidence[n]
            assert all(Q[(n,u)]==binary_windows[(L,n,u)] for u in range(2,L+1))
            assert D[n]==r['Delta'] and B_from_Q(L,n,lambda _L,u,_n:Q[(n,u)])[0]==r['B_real']
            for field in ('all_optimal_subsets','anchored_optima','translation_reflection_classes'):
                assert e[field]==r[field],(L,n,field,e[field],r[field])
            assert e['witness']==r['witness']
            brute_pairs+=1;optimum_count_pairs+=1
        print('Enumerated interval',L,flush=True)
    report=dict(status='PASS',paper='DAM',finite_witnesses=252,boundary_tight=245,binary_tight=247,
        complete_enumerated_pairs=brute_pairs,optimum_counts_and_classes_checked=optimum_count_pairs,
        moment_hierarchy_vectors=5,joint_integer_vectors=714,candidate_reduction_fixtures=verify_reduction(),
        canonical_table_exports=regenerate(),full_map_rerun=full,seconds=round(time.perf_counter()-started,3),
        prespecified_grid_pairs=462,strict_boundary_over_integer_improvements=189,principal_medium_witnesses=9,
        complete_walk_windows=walk_windows, strengthened_hierarchy_vectors=5,modular_construction_runs=len(modular),
        direct_coverage_split=splits['direct'],transferred_coverage_split=splits['transferred'],trivial_rows=21,
        binary_window_minima_checked=sum(L-1 for L in range(4,25 if full else 13) for n in range(2,L)),
        scope='Finite computation checks; not a replacement for human review of the universal proofs.')
    out=HERE/'recomputed';out.mkdir(exist_ok=True)
    (out/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--full',action='store_true');a=p.parse_args()
    print(json.dumps(check(a.full),indent=2))
