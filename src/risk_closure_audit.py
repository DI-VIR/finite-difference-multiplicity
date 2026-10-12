"""Symbolic identities and small MILP checks, not independent human proof review."""
from pathlib import Path
import json
import sympy as sp
from binary_window_milp import validate

ROOT = Path(__file__).resolve().parents[1]


def audit():
    L, n, u, d = sp.symbols('L n u d', integer=True, positive=True)
    k = sp.symbols('k', integer=True, positive=True)
    t = sp.symbols('t', positive=True)
    identities = []

    def check(name, left, right):
        assert sp.simplify(left - right) == 0, name
        identities.append(name)

    check('unordered-pair weighted expectation',
          2 * sp.summation((u - d) * (L - d), (d, 1, u - 1)),
          u * (u - 1) * (L - (u + 1) / 3))
    raw_cs = ((n * u)**2 / (L + u - 1) - n * u) / (u * (u - 1))
    check('half-density real expression', raw_cs.subs({n: L / 2, u: t + 1}),
          -L / (4 * t) + L * (L - 1) / (4 * (L + t)))
    check('integer correction split', (L + u - 1) / (4 * u * (u - 1)),
          L / (4 * u * (u - 1)) + 1 / (4 * u))
    check('small-window upper expression',
          L / 4 - sp.Rational(1, 4) - L / (4 * (u - 1))
          + L / (4 * u * (u - 1)) + 1 / (4 * u),
          L / 4 - sp.Rational(1, 4) - (L - 1) / (4 * u))
    layer = k - 1
    mass = 2 * k**2 * (2 * k)
    eta = sp.expand(mass - layer * (4 * k**2 + 2 * k - layer))
    check('u=2k residual', eta, 3 * k**2 + 1)
    q = ((layer - 1) * layer * (2 * layer - 1) / 3
         + (4 * k**2 + 2 * k - 2 * layer + 1) * layer**2
         + (2 * layer + 1) * eta)
    check('u=2k capacity polynomial', q,
          4 * k**4 - sp.Rational(4, 3) * k**3 + k**2 + k / 3)
    check('u=2k raw boundary expression', (q - mass) / (2 * k * (2 * k - 1)),
          k**2 - sp.Rational(5, 6) * k - sp.Rational(1, 6))
    check('u=2k+1 real expression', raw_cs.subs({L: 4 * k**2, n: 2 * k**2, u: 2 * k + 1}),
          k**2 - k)
    for residue in range(6):
        j = sp.symbols('j', integer=True)
        kval = 6 * j + residue
        check('ceiling/floor residue %d' % residue,
              kval - sp.floor((5 * kval + 1) / 6) - 2,
              sp.ceiling((kval - 1) / 6) - 2)
    validation = validate()
    saved = json.loads((ROOT / 'data/binary_window_milp_validation.json').read_text())
    assert [(r['L'], r['n'], r['u'], r['Q_real']) for r in validation['records']] == [
        (r['L'], r['n'], r['u'], r['Q_real']) for r in saved['records']]
    return dict(status='PASS', symbolic_identity_checks=identities,
                weighted_binary_milp=validation,
                proof_status='Auxiliary symbolic and finite checks only; human audit pending.')


if __name__ == '__main__':
    report = audit()
    (ROOT / 'tests/risk_closure_audit.json').write_text(
        json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
