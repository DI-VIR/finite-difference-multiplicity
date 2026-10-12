"""Reject incomplete claim coverage and invalid Fourier input domains."""
from pathlib import Path
import copy
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fourier_verify import verify_packet, run

trees=json.loads((ROOT/'data/fourier_exact_trees.json').read_text())
comparisons=json.loads((ROOT/'data/fourier_comparison.json').read_text())
profiles={(p['L'],p['n']):p['per_u'] for p in json.loads((ROOT/'data/seven_case_profiles.json').read_text())}


def rejected(t, c, p):
    try:
        verify_packet(t,c,p)
    except AssertionError:
        return
    raise AssertionError('A damaged claim packet was accepted.')


verify_packet(trees,comparisons,profiles)
rejected(trees[:-1],comparisons,profiles)
rejected(trees+[trees[0]],comparisons,profiles)
rejected(trees,comparisons[:-1],profiles)
bad=copy.deepcopy(trees)
bad[0]['cap']+=1
rejected(bad,comparisons,profiles)
bad=copy.deepcopy(trees)
bad[0]['rational_cos_grid'][0]='-129/128'
rejected(bad,comparisons,profiles)
bad=copy.deepcopy(profiles)
bad[(17,6)]=bad[(17,6)][:-1]
rejected(trees,comparisons,bad)
bad=copy.deepcopy(profiles)
bad[(17,6)].append(bad[(17,6)][0])
rejected(trees,comparisons,bad)
assert run()['cases']==5
print('PASS: seven packet mutations rejected; all five exact exclusions verified.')
