"""Check every saved exact pair and the full 462-pair integer-bound grid."""
from _bootstrap import load,report
from extremal import bounds

def check():
    records=load('exact_map_nontrivial')
    assert len(records)==252
    for row in records:
        b=bounds(row['L'],row['n'])
        assert b['B_bd']==row['B_bd']<=row['B_real']<=row['Delta']
    grid=load('bound_landscape')
    assert len(grid)==462
    for row in grid:
        b=bounds(row['L'],row['n'])
        assert all(b[k]==row[k] for k in ('B_ph','B_CS','B_int','B_bd','u_bd'))
        assert row['bd_minus_int']==b['B_bd']-b['B_int']
    assert sum(row['bd_minus_int']>0 for row in grid)==189
    report(exact_pairs=252,grid_pairs=462,strict_boundary_improvements=189)

if __name__=='__main__':check()
