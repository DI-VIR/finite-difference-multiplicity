"""Export additional comparisons and audited historical witness provenance."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import argparse
import csv
import json

from extremal import bounds
from regenerate_tables import read_tables, write_tables
ROOT=Path(__file__).resolve().parents[1]


def load(name): return json.loads((ROOT/'data'/(name+'.json')).read_text())


def save(name,data):
    (ROOT/'data'/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n')
    if data and isinstance(data,list) and isinstance(data[0],dict):
        keys=list(dict.fromkeys(k for r in data for k in r))
        with (ROOT/'data'/(name+'.csv')).open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=keys,lineterminator='\n');w.writeheader()
            w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in data)


def xu_bound(L,n):
    # Printed Theorem 3.1 and Proposition 2.2 of arXiv:1706.05969v5.
    candidates=[]
    for D in range(2,L+1):
        fourier=Fraction(n*n,L+D)-Fraction(n+1,D)
        inclusion=Fraction(2*(D*n-L),D*(D-1))-Fraction(2,D)
        for label,value in [('Theorem 3.1',fourier),('Proposition 2.2',inclusion)]:
            candidates.append((-((-value.numerator)//value.denominator),D,label))
    value,D,label=max(candidates,key=lambda r:(r[0],-r[1]))
    return dict(B_Xu=max(0,value),D_Xu=D,Xu_source=label)


def import_history(workspace):
    paths=[workspace/'ablation_results/20260920_comprehensive/runs.json',
           workspace/'ablation_results/20260923_milp_optimum_certificates/runs.json',
           workspace/'ablation_results/20260923_milp_optimum_certificates/extended_upper_attempts.json']
    sources=[(p,json.loads(p.read_text()),sha256(p.read_bytes()).hexdigest()) for p in paths]
    protocol=json.loads((paths[0].parent/'protocol.json').read_text())
    records=[]
    for medium in load('verified_medium'):
        L,n=medium['L'],medium['n']
        taps=list(map(int,medium['witness'].split()))
        matches=[]
        for path,rows,digest in sources:
            for index,row in enumerate(rows):
                original=row.get('taps') or row.get('verified_taps')
                if (row.get('L'),row.get('n'))==(L,n) and original is not None and sorted(original)==taps:
                    matches.append((path,index,row,digest))
        assert matches,(L,n)
        path,index,row,digest=matches[0]
        construction=path==paths[0]
        method=row.get('method',row.get('solver','SciPy/HiGHS feasibility'))
        record=dict(L=L,n=n,value=medium['best_feasible'],witness=taps,
            source=path.relative_to(workspace).as_posix(),source_record_index=index,
            original_file_sha256=digest,source_kind='Historical construction' if construction else 'Historical solver',
            method=method,method_label=row.get('method_label',method),
            seed=row.get('seed'),proposal_budget=row.get('eval_budget'),
            optimization_seed=row['seed']+100000+7919*protocol['methods'].index(method) if construction and method!='random' else None,
            witness_sha256=sha256(json.dumps(taps,separators=(',',':')).encode()).hexdigest(),
            solver_witness_is_same=any(r['L']==L and r['n']==n and r['witness']==taps for r in load('milp_complete')),
            historical_not_new_modular_arm=True)
        records.append(record)
    save('medium_witness_provenance',records)
    save('historical_construction_protocol',dict(
        source='tools/comprehensive_tap_experiments.py',
        source_sha256=sha256((workspace/'tools/comprehensive_tap_experiments.py').read_bytes()).hexdigest(),
        original_run_protocol={k:protocol[k] for k in ('parameters','seeds','evals_per_optimizer','methods','random_seed_start')},
        LB_SOM_objective='1000 sum_d max(r_d-B_int,0)^2 + 10 max_d r_d',
        note='Only geometry witness fields are extracted; no cryptographic endpoints are reproduced or used as evidence in this paper.'))


def run(workspace=None):
    if workspace: import_history(workspace)
    domain={(r['L'],r['n']) for r in load('literature_grid_comparison')}
    domain.update((r['L'],r['n']) for r in load('exact_map_nontrivial'))
    domain.update((r['L'],r['n']) for r in load('verified_medium'))
    comparisons=[]
    for L,n in sorted(domain):
        b=bounds(L,n)
        comparisons.append(dict(L=L,n=n,**xu_bound(L,n),B_CS=b['B_CS'],B_int=b['B_int'],B_bd=b['B_bd']))
    save('xu_comparison',comparisons)
    cat=read_tables()
    def add(key,caption,headers,rows):
        assert all(len(r)==len(headers) for r in rows)
        cat[key]=dict(caption=caption,headers=headers,rows=rows)
    checked=load('fourier_exact_trees')
    eigen=load('fourier_comparison')
    exact={(r['L'],r['n']):r['Delta'] for r in load('exact_map_nontrivial')}
    add('S16','Five exact exclusions after Fourier strengthening. Rational negative values invalidate the old vectors. Complete branch trees with exact rational inequality certificates exclude every integer vector at the old cap, not only those old vectors. Finite L-by-L Toeplitz eigenvalues are numerical diagnostics; the tree exclusions are exact.',
        ['(L,n)','Old cap',r'$\Delta$','Old-vector F value','Toeplitz min eigenvalue','Tree leaves','Fourier cap optimum'],
        [[f"({r['L']},{r['n']})",r['cap'],exact[(r['L'],r['n'])],f"{e['negative_value']:.4f}",f"{e['finite_Toeplitz_min_eigenvalue']:.4f}",sum(len(b['leaves']) for b in r['branches']),exact[(r['L'],r['n'])]] for r,e in zip(checked,eigen)])
    medium=load('verified_medium')
    cells={(r['L'],r['n']):r for r in comparisons}
    add('S17','Same-model finite comparison with Xu [22]. B_Xu maximizes the ceilings of the printed Theorem 3.1 and Proposition 2.2 over 2<=D<=L; no asymptotic remainder or unproved constant is inserted. These are lower certificates, not heuristic upper witnesses.',
        ['(L,n)',r'$B_{\mathrm{Xu}}$','B_CS','B_int','B_bd','Verified upper'],
        [[f"({m['L']},{m['n']})",cells[(m['L'],m['n'])]['B_Xu'],m['B_CS'],m['B_int'],m['B_bd'],m['best_feasible']] for m in medium])
    summaries=load('seed_budget_summary')
    rows=[]
    for L,n in [(m['L'],m['n']) for m in medium]:
        for family in ('Paley','Singer','Ruzsa-Spence','RDS','Random','Random-SA'):
            row=[f'({L},{n})',family]
            for budget in (0,3000,30000):
                s=next(r for r in summaries if (r['L'],r['n'],r['family'],r['budget'])==(L,n,family,budget))
                row.extend([f"{s['mean']:.2f} +/- {s['SD']:.2f}",s['best']])
            rows.append(row)
    add('S18','Initialization and budget ablation. Each cell has 30 streams. Columns show mean plus sample SD and best multiplicity. Zero budget is seed only, not a separate random population. Random-SA retains the best incumbent under the disclosed annealing schedule; all arms optimize direct multiplicity.',
        ['(L,n)','Arm','0 mean +/- SD','0 best','3k mean +/- SD','3k best','30k mean +/- SD','30k best'],rows)
    add('S19','Exploratory paired contrasts at 30000 proposals, conditional on the declared deterministic algebraic seeds. Difference is algebraic minus random hill climbing, so negative favors the algebraic arm. Pointwise percentile bootstrap intervals resample 30 shared seed blocks 10000 times; they are not familywise-adjusted superiority tests. Both budgets and full raw outputs are archived.',
        ['(L,n)','Algebraic arm','Mean difference','95% lower','95% upper'],
        [[f"({r['L']},{r['n']})",r['family'],f"{r['algebraic_minus_random_mean']:.3f}",f"{r['paired_bootstrap_95_low']:.3f}",f"{r['paired_bootstrap_95_high']:.3f}"] for r in load('seed_budget_contrasts') if r['budget']==30000])
    windows=load('window_scalability')
    for r in windows:
        from binary_window_milp import square_sum
        assert square_sum(r['witness'],r['u'])==r['Q_upper']
    add('S20','Targeted binary-window MILP beyond L=24. Each run has a 60-second limit. Q_bd is an exact analytic lower bound; the dual column is the numerical solver lower bound. The upper column is checked by exact unordered-pair counting. Optimal refers to solver status, not an exported exact dual certificate. A time-limit bracket does not compute B_real over all windows.',
        ['(L,n,u)','Q_bd','Solver dual','Checked upper','Outcome','Seconds'],
        [[f"({r['L']},{r['n']},{r['u']})",__import__('extremal').min_boundary_square(r['L'],r['n'],r['u'])[0],f"{r['Q_lower_numerical']:.3f}",r['Q_upper'],'Optimal' if r['status']==0 else 'Time limit',f"{r['seconds']:.2f}"] for r in windows])
    provenance=load('medium_witness_provenance')
    add('S21','Exact source matching of the nine historical upper witnesses. These are distinct from Table 9 modular-search outputs. Record indices are zero based. Complete coordinates, original source-file hashes and witness hashes are archived. Separate HiGHS witnesses on overlapping pairs need not have the same coordinates.',
        ['(L,n)','Value','Source kind','Method','Seed','Proposals','Source row'],
        [[f"({r['L']},{r['n']})",r['value'],r['source_kind'],r['method_label'],r['seed'] if r['seed'] is not None else 'Not recorded',r['proposal_budget'] if r['proposal_budget'] is not None else 'Not recorded',r['source_record_index']] for r in provenance])
    for row in cat['S15']['rows']:
        if row[0]=='B_real benchmark': row[2]=r'Universal validity; full map $L\leq24$; six targeted MILP runs $28\leq L\leq80$'
        if row[0]=='Fixed-cap NP-completeness': row[0]='Fixed-cap strong NP-completeness'
    cat['7']['rows'][-1][-1]='Strongly NP-complete for each fixed g (Theorem 10)'
    cat['S7']['rows']=[r for r in cat['S7']['rows'] if r[0]!='Xu [22]']
    cat['S7']['rows'].append(['Xu [22]','Popular positive differences','Finite N, density and threshold D',
        'Finite Fourier and inclusion-exclusion lower bounds','Yes, with D=L','Direct antecedent; finite comparison in Section S17'])
    write_tables(cat)
    print('Built',len(cat),'tables; added S16 through S21.')
    print('Historical sources:',[(r['L'],r['n'],r['method'],r['seed']) for r in provenance])
    grid=[r for r in comparisons if (r['L'],r['n']) in {(g['L'],g['n']) for g in load('literature_grid_comparison')}]
    print('Xu vs B_CS:',Counter('lower' if r['B_Xu']<r['B_CS'] else 'equal' if r['B_Xu']==r['B_CS'] else 'higher' for r in grid))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--import-history',type=Path)
    args=parser.parse_args();run(args.import_history)
