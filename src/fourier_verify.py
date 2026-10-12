"""Solver-independent exact checking of stored Fourier exclusion trees."""
from fractions import Fraction as F
from pathlib import Path
import json
from fourier_comparison import coefficients, exact_value

ROOT=Path(__file__).resolve().parents[1]
EXPECTED_CASES={(17,6):1,(18,8):2,(20,11):4,(20,12):5,(23,9):2}


def verify_packet(trees, comparisons, profiles):
    """Bind verification to the complete five-case claim, not a file subset."""
    keys=[(t['L'],t['n']) for t in trees]
    assert len(keys)==len(set(keys))==len(EXPECTED_CASES)
    assert set(keys)==set(EXPECTED_CASES)
    comparison_keys=[(p['L'],p['n']) for p in comparisons]
    assert len(comparison_keys)==len(set(comparison_keys))==len(EXPECTED_CASES)
    assert set(comparison_keys)==set(EXPECTED_CASES)
    for tree in trees:
        L,n=tree['L'],tree['n']
        assert tree['cap']==EXPECTED_CASES[(L,n)]
        cuts=[F(x) for x in tree['rational_cos_grid']]
        assert cuts and len(cuts)==len(set(cuts))
        assert all(F(-1)<=x<=F(1) for x in cuts)
        rows=profiles[(L,n)]
        assert sorted(p['u'] for p in rows)==list(range(2,L+1))
        assert all(isinstance(p['Q_real'],int) and p['Q_real']>=n*p['u'] for p in rows)
    for comparison in comparisons:
        assert F(-1)<=F(comparison['cos_theta'])<=F(1)


def verify_leaf(L,n,g,mass,profile,cuts,leaf):
    degree=L-1
    lo,hi=leaf['lower'],leaf['upper']
    assert len(lo)==len(hi)==degree
    assert all(isinstance(v,int) for v in lo+hi)
    assert all(0<=a<=b<=min(g,n,L-d) for d,(a,b) in enumerate(zip(lo,hi),1))
    base=[[F(1)]*degree]
    lower=[n*(n-1)//2]; upper=lower.copy()
    for p in profile:
        base.append([F(2*max(0,p['u']-d)) for d in range(1,L)])
        lower.append(p['Q_real']-n*p['u']); upper.append(None)
    for j in range(degree):
        base.append([F(int(i==j)) for i in range(degree)])
        lower.append(lo[j]); upper.append(hi[j])
    A,b=[],[]
    for row,low,high in zip(base,lower,upper):
        A.append([-v for v in row]); b.append(F(-low))
        if high is not None: A.append(row); b.append(F(high))
    odd=[F(int(d%2)) for d in range(1,L)]
    A.extend([odd,[-v for v in odd]]); b.extend([F(mass),F(-mass)])
    for row in cuts: A.append([-v for v in row]); b.append(F(n))
    for i,limit in enumerate(hi):
        row=[F(int(j==i)) for j in range(degree)]
        A.extend([row,[-v for v in row]]); b.extend([F(limit),F(0)])
    cert=leaf['certificate']
    active=[(int(i),F(w)) for i,w in cert['weights']]
    assert len({i for i,_ in active})==len(active) and all(0<=i<len(A) and w>0 for i,w in active)
    total=sum(w*b[i] for i,w in active)
    residual=[sum(w*A[i][j] for i,w in active) for j in range(degree)]
    floor=sum(min(v*limit,F(0)) for v,limit in zip(residual,hi))
    assert total<floor and str(total)==cert['combined_rhs'] and str(floor)==cert['residual_lower_bound']


def check_cover(leaves,lo,hi,prefix=()):
    here=[p for p in leaves if tuple(tuple(v) for v in p['path'][:len(prefix)])==prefix]
    assert here
    terminal=[p for p in here if len(p['path'])==len(prefix)]
    if terminal:
        assert len(here)==len(terminal)==1
        assert terminal[0]['lower']==lo and terminal[0]['upper']==hi
        return
    next_steps={tuple(p['path'][len(prefix)]) for p in here}
    assert len(next_steps)==2
    left=next(s for s in next_steps if s[1]=='<=')
    right=next(s for s in next_steps if s[1]=='>=')
    j,_,k=left
    assert right==(j,'>=',k+1) and lo[j]<=k<hi[j]
    lhi=hi.copy(); lhi[j]=k
    rlo=lo.copy(); rlo[j]=k+1
    check_cover(here,lo.copy(),lhi,prefix+(left,))
    check_cover(here,rlo,hi.copy(),prefix+(right,))


def run():
    comparisons=json.loads((ROOT/'data/fourier_comparison.json').read_text())
    for p in comparisons:
        r=next(v for v in json.loads((ROOT/'data/residual_diagnostics.json').read_text()) if (v['L'],v['n'])==(p['L'],p['n']))['strengthened_witness']
        assert exact_value(p['n'],r,F(p['cos_theta']))==F(p['exact_negative_value'])<0
    profiles={(p['L'],p['n']):p['per_u'] for p in json.loads((ROOT/'data/seven_case_profiles.json').read_text())}
    trees=json.loads((ROOT/'data/fourier_exact_trees.json').read_text())
    verify_packet(trees,comparisons,profiles)
    counts=[]
    for tree in trees:
        L,n,g=tree['L'],tree['n'],tree['cap']
        assert tree['outcome']=='EXACT_RATIONAL_INTEGER_EXCLUSION'
        expected=sorted({e*(n-e) for e in range(max(0,n-L//2),min(n,(L+1)//2)+1)})
        assert [p['odd_mass'] for p in tree['branches']]==expected
        cuts=[coefficients(F(x),L-1) for x in tree['rational_cos_grid']]
        count=0
        for branch in tree['branches']:
            assert branch['complete'] and branch['remaining_nodes']==0
            check_cover(branch['leaves'],[0]*(L-1),[min(g,n,L-d) for d in range(1,L)])
            for leaf in branch['leaves']:
                verify_leaf(L,n,g,branch['odd_mass'],profiles[(L,n)],cuts,leaf)
                count+=1
        counts.append(dict(L=L,n=n,cap=g,exact_leaves=count))
    report=dict(status='PASS',cases=len(counts),records=counts,
        method='Exact Fraction inequalities and complete integer split coverage; no optimization call in this verification entrypoint.')
    (ROOT/'tests/fourier_verify.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return report


if __name__=='__main__': run()
