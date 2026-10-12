"""Verify all five residual moment vectors and 714 joint certificates."""
from _bootstrap import load,report
from realizability_bound import Q_bd

def check():
    residual=[row for row in load('realizability_audit') if row.get('hierarchy_exhausted')]
    assert len(residual)==5
    for row in residual:
        vec=row['integer_r_at_B_bd']
        assert len(vec)==row['L']-1 and min(vec)>=0
        assert sum(vec)==row['n']*(row['n']-1)//2
        assert max(vec)==row['B_bd']<row['Delta']
        for p in row['per_u']:
            assert 2*sum((p['u']-d)*vec[d-1] for d in range(1,p['u']))>=p['Q_real']-row['n']*p['u']
    joint=load('joint_bound_audit')
    assert len(joint)==714
    for row in joint:
        L,n=row['L'],row['n'];v=row['integer_witness']
        assert min(v)>=0 and max(v)<=row['max_bd_ph'] and sum(v)==n*(n-1)//2
        assert all(2*sum((u-d)*v[d-1] for d in range(1,u))>=Q_bd(L,u,n)-n*u for u in range(2,L+1))
    report(residual_vectors=5,joint_vectors=714,scope='Finite certificate checks, not a universal proof')

if __name__=='__main__':check()
