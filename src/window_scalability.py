"""Targeted L>24 weighted binary MILP runs with explicit unresolved outcomes."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix
from binary_window_milp import square_sum

ROOT=Path(__file__).resolve().parents[1]
INSTANCES=((28,12,7),(32,15,8),(40,19,10),(48,22,12),(64,29,16),(80,37,20))


def solve(L,n,u,limit=60):
    products=[(i,d) for d in range(1,u) for i in range(L-d)]
    size=L+len(products)
    objective=np.zeros(size)
    rows,cols,values=[],[],[]
    lo,hi=[float(n)],[float(n)]
    for i in range(L): rows.append(0); cols.append(i); values.append(1.0)
    for j,(i,d) in enumerate(products,L):
        objective[j]=u-d
        for entries,low,high in ((((j,1),(i,-1)),-np.inf,0),
                (((j,1),(i+d,-1)),-np.inf,0),(((j,1),(i,-1),(i+d,-1)),-1,np.inf)):
            r=len(lo)
            for c,v in entries: rows.append(r); cols.append(c); values.append(v)
            lo.append(low); hi.append(high)
    matrix=coo_matrix((values,(rows,cols)),shape=(len(lo),size)).tocsc()
    tick=time.perf_counter()
    result=milp(objective,integrality=np.r_[np.ones(L),np.zeros(len(products))],
        bounds=Bounds(np.zeros(size),np.ones(size)),constraints=LinearConstraint(matrix,lo,hi),
        options={'time_limit':float(limit),'mip_rel_gap':0.0})
    record=dict(L=L,n=n,u=u,time_limit=limit,seconds=round(time.perf_counter()-tick,6),
        status=int(result.status),message=result.message,variables=size,constraints=len(lo),
        nodes=int(getattr(result,'mip_node_count',0) or 0),
        Q_lower_numerical=n*u+2*float(result.mip_dual_bound) if getattr(result,'mip_dual_bound',None) is not None else None,
        Q_upper=None,witness=None,interpretation='One window only; numerical solver certificate, independently checked incumbent.')
    if result.x is not None:
        marks=np.flatnonzero(np.rint(result.x[:L]).astype(int)).tolist()
        assert len(marks)==n and np.allclose(result.x[:L],np.rint(result.x[:L]),atol=1e-6)
        q=square_sum(marks,u)
        assert abs(result.fun-(q-n*u)/2)<1e-5
        record.update(Q_upper=q,witness=marks)
    record['outcome']='SOLVER_OPTIMAL' if result.status==0 else 'TIME_LIMIT_BRACKET' if result.status==1 else 'UNRESOLVED'
    return record


if __name__=='__main__':
    records=[]
    for args in INSTANCES:
        r=solve(*args)
        records.append(r)
        (ROOT/'data/window_scalability.json').write_text(json.dumps(records,indent=2)+'\n')
        print(args,r['outcome'],r['Q_lower_numerical'],r['Q_upper'],r['seconds'],flush=True)
