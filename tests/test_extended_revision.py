"""New science records and rejection of damaged exact certificate trees."""
from pathlib import Path
import copy
import json
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from extended_audit import audit
from fourier_verify import check_cover, verify_leaf
from fourier_comparison import coefficients
from fractions import Fraction

audit()
trees=json.loads((ROOT/'data/fourier_exact_trees.json').read_text())
tree=trees[1]
branch=next(b for b in tree['branches'] if len(b['leaves'])>1)
try:
    check_cover(branch['leaves'][:-1],[0]*(tree['L']-1),
                [min(tree['cap'],tree['n'],tree['L']-d) for d in range(1,tree['L'])])
except AssertionError:
    pass
else:
    raise AssertionError('A missing exclusion leaf was not rejected.')
damaged=copy.deepcopy(tree['branches'][0]['leaves'][0])
damaged['certificate']['weights'][0][1]='-1'
profiles=json.loads((ROOT/'data/seven_case_profiles.json').read_text())
profile=next(p['per_u'] for p in profiles if (p['L'],p['n'])==(tree['L'],tree['n']))
cuts=[coefficients(Fraction(x),tree['L']-1) for x in tree['rational_cos_grid']]
try:
    verify_leaf(tree['L'],tree['n'],tree['cap'],tree['branches'][0]['odd_mass'],profile,cuts,damaged)
except AssertionError:
    pass
else:
    raise AssertionError('A negative certificate weight was not rejected.')
print('PASS: independent record audit and damaged tree/weight mutation rejection.')
