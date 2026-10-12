"""Necessary tail-capacity and parity conditions, tested without claiming sufficiency."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

HERE = Path(__file__).resolve().parents[1]


def run():
    records = json.loads((HERE/'data/realizability_audit.json').read_text())
    output = []
    for r in records:
        if not r.get('hierarchy_exhausted'):
            continue
        L,n,g = r['L'],r['n'],r['B_bd']
        v = r['integer_r_at_B_bd']
        odd_mass = sum(v[::2])
        possible = sorted({e*(n-e) for e in range(n+1)})
        tail_failures = [d for d in range(1,L) if v[d-1] > min(g,n,L-d)]
        A = [np.ones(L-1)]
        lows = [n*(n-1)/2]
        highs = lows.copy()
        for p in r['per_u']:
            A.append(np.array([2*max(0,p['u']-d) for d in range(1,L)],dtype=float))
            lows.append(p['Q_real']-n*p['u'])
            highs.append(np.inf)
        parity = np.zeros(L-1);parity[::2] = 1
        matrix = np.vstack(A+[parity])
        ub = np.array([min(g,n,L-d) for d in range(1,L)])
        statuses, witness = [], None
        for mass in possible:
            result = milp(c=np.zeros(L-1),integrality=np.ones(L-1),bounds=Bounds(0,ub),
                constraints=LinearConstraint(matrix,np.array(lows+[mass]),np.array(highs+[mass])),
                options={'time_limit':30})
            statuses.append(dict(odd_mass=mass,status=int(result.status),message=result.message))
            if result.status == 0:
                w = np.rint(result.x).astype(int).tolist()
                assert sum(w[::2]) == mass and sum(w) == n*(n-1)//2
                assert all(0 <= a <= b for a,b in zip(w,ub))
                assert all(2*sum((p['u']-d)*w[d-1] for d in range(1,p['u'])) >= p['Q_real']-n*p['u'] for p in r['per_u'])
                witness = w
                break
        output.append(dict(L=L,n=n,cap=g,old_odd_mass=odd_mass,possible_odd_mass=possible,
            old_parity_feasible=odd_mass in possible,old_tail_violations=tail_failures,
            strengthened_status='FEASIBLE' if witness else 'INFEASIBLE' if all(x['status']==2 for x in statuses) else 'UNKNOWN',
            strengthened_witness=witness,solver_cases=statuses))
    (HERE/'data/residual_diagnostics.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(output,indent=2))


if __name__ == '__main__':
    run()
