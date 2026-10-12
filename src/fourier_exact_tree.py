"""Exact validation of a finite Fourier-cut branch tree, when it terminates."""
from fractions import Fraction
from pathlib import Path
import json
import time

import numpy as np
from scipy.optimize import linprog
from fourier_comparison import coefficients, rational_exclusion

ROOT=Path(__file__).resolve().parents[1]


def prove(case, profiles, seconds=180, max_nodes=10000):
    L,n,g=case['L'],case['n'],case['cap']
    degree=L-1
    cuts=[Fraction(i,128) for i in range(-128,129)]
    comparisons=json.loads((ROOT/'data/fourier_comparison.json').read_text())
    comparison=next(p for p in comparisons if (p['L'],p['n'])==(L,n))
    cuts.extend(sorted({Fraction(a['negative_cut']) for b in comparison['new_relaxation']['parity_branches']
                        for a in b['attempts'] if 'negative_cut' in a}-set(cuts)))
    matrix=np.asarray([[float(v) for v in coefficients(x,degree)] for x in cuts])
    moments=np.asarray([[2*max(0,p['u']-d) for d in range(1,L)] for p in profiles])
    moment_rhs=np.asarray([p['Q_real']-n*p['u'] for p in profiles])
    A=-np.vstack([moments,matrix])
    b=-np.r_[moment_rhs,np.full(len(cuts),-n)]
    odd=np.asarray([int(d%2) for d in range(1,L)])
    E=np.asarray([np.ones(degree),odd])
    masses=sorted({e*(n-e) for e in range(max(0,n-L//2),min(n,(L+1)//2)+1)})
    started=time.perf_counter()
    nodes=0
    branches=[]
    for mass in masses:
        leaves=[]
        stack=[([0]*degree,[min(g,n,L-d) for d in range(1,L)],[])]
        complete=True
        while stack:
            if time.perf_counter()-started>seconds or nodes>=max_nodes:
                complete=False; break
            lo,hi,path=stack.pop()
            nodes+=1
            result=linprog(np.zeros(degree),A_ub=A,b_ub=b,A_eq=E,
                b_eq=[n*(n-1)/2,mass],bounds=list(zip(lo,hi)),method='highs')
            if result.status==2:
                base=[np.ones(degree),*moments]
                lower=[n*(n-1)//2,*moment_rhs]
                upper=[n*(n-1)//2,*([np.inf]*len(moments))]
                for j in range(degree):
                    row=np.asarray([int(i==j) for i in range(degree)])
                    base.append(row); lower.append(lo[j]); upper.append(hi[j])
                cert=rational_exclusion(base,lower,upper,cuts,n,hi,mass)
                if cert['outcome']!='EXACT_RATIONAL_EXCLUSION':
                    complete=False; break
                cert.pop('cuts')
                leaves.append(dict(path=path,lower=lo,upper=hi,certificate=cert))
                continue
            if result.status!=0:
                complete=False; break
            differences=abs(result.x-np.rint(result.x))
            j=int(np.argmax(differences))
            if differences[j]<1e-7:
                complete=False; break
            floor=int(np.floor(result.x[j]))
            left=hi.copy(); left[j]=min(left[j],floor)
            right=lo.copy(); right[j]=max(right[j],floor+1)
            if lo[j]<=left[j]: stack.append((lo.copy(),left,path+[[j,'<=',floor]]))
            if right[j]<=hi[j]: stack.append((right,hi.copy(),path+[[j,'>=',floor+1]]))
        branches.append(dict(odd_mass=mass,complete=complete,leaves=leaves,
            remaining_nodes=len(stack)))
    return dict(L=L,n=n,cap=g,nodes=nodes,seconds=round(time.perf_counter()-started,6),
        outcome='EXACT_RATIONAL_INTEGER_EXCLUSION' if all(p['complete'] for p in branches) else 'INCOMPLETE_EXACT_TREE',
        rational_cos_grid=[str(x) for x in cuts],branches=branches)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--only',type=int,default=0)
    args=parser.parse_args()
    cases=json.loads((ROOT/'data/residual_diagnostics.json').read_text())
    profiles={(p['L'],p['n']):p['per_u'] for p in json.loads((ROOT/'data/seven_case_profiles.json').read_text())}
    output=ROOT/'data/fourier_exact_trees.json'
    records=json.loads(output.read_text()) if args.only else []
    for case in cases:
        if args.only and case['L']!=args.only: continue
        report=prove(case,profiles[(case['L'],case['n'])])
        records=[p for p in records if (p['L'],p['n'])!=(case['L'],case['n'])]+[report]
        records.sort(key=lambda p:(p['L'],p['n']))
        (ROOT/'data/fourier_exact_trees.json').write_text(json.dumps(records,indent=2)+'\n')
        print(report['L'],report['n'],report['outcome'],report['nodes'],report['seconds'],flush=True)
