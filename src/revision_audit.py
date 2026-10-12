"""Additional exact checks for the submission revision; no stochastic proof claims."""
from collections import Counter
from pathlib import Path
import csv
import json
import time

from extremal import bounds, min_boundary_square
from realizability_bound import Q_walk

HERE = Path(__file__).resolve().parents[1]


def integer_walk(L, n, u):
    """Independent sparse integer DP, also returning an attaining profile."""
    caps = [min(L - 1, x + u - 1) - max(0, x) + 1 for x in range(1 - u, L)]
    M = n * u
    states = {(a, a): (a * a, (a,)) for a in (0, 1)}
    for i, cap in enumerate(caps[1:], 1):
        nxt = {}
        for (mass, last), (cost, path) in states.items():
            for a in range(max(0, last - 1), min(cap, last + 1) + 1):
                if mass + a > M or (i == len(caps) - 1 and a > 1):
                    continue
                key, value = (mass + a, a), cost + a * a
                if key not in nxt or value < nxt[key][0]:
                    nxt[key] = (value, path + (a,))
        states = nxt
    cost, profile = min(v for (mass, last), v in states.items() if mass == M and last <= 1)
    assert sum(profile) == M and sum(a * a for a in profile) == cost
    assert all(0 <= a <= h for a, h in zip(profile, caps))
    assert all(abs(a - b) <= 1 for a, b in zip(profile, profile[1:]))
    return cost, list(profile)


def run():
    started = time.perf_counter()
    audit = json.loads((HERE / 'data/realizability_audit.json').read_text())
    rows, profiles, signatures = [], [], []
    for record in audit:
        if record['B_bd'] == record['Delta']:
            continue
        L, n = record['L'], record['n']
        per_u = []
        for old in record['per_u']:
            u = old['u']
            qwalk, path = integer_walk(L, n, u)
            assert qwalk == Q_walk(L, u, n)
            qbd = min_boundary_square(L, n, u)[0]
            assert qbd == old['Q_bd'] <= qwalk <= old['Q_real']
            per_u.append(dict(u=u, Q_bd=qbd, Q_walk=qwalk, Q_real=old['Q_real'], walk_profile=path))
        result = dict(L=L, n=n, Delta=record['Delta'])
        for layer in ('bd', 'walk', 'real'):
            values = [(-(-(p['Q_' + layer] - n * p['u']) // (p['u'] * (p['u'] - 1))), p['u']) for p in per_u]
            best = max(v for v, _ in values)
            args = [u for v, u in values if v == best]
            result['B_' + layer] = best
            result['u_' + layer] = min(args)
            result['argmax_' + layer] = args
        ustar = result['u_real']
        point = next(p for p in per_u if p['u'] == ustar)
        result.update({key: point[key] for key in ('Q_bd', 'Q_walk', 'Q_real')})
        result['Q_reference_window'] = ustar
        result['gap_closed'] = result['B_real'] == result['Delta']
        rows.append(result)
        profiles.append(dict(L=L, n=n, per_u=per_u))
        record.update(B_walk=result['B_walk'], u_walk=result['u_walk'])
        for p in record['per_u']:
            p['Q_walk'] = next(q['Q_walk'] for q in per_u if q['u'] == p['u'])
        vec = record.get('integer_r_at_B_bd')
        if vec is not None and record.get('hierarchy_exhausted'):
            slacks = [2 * sum((p['u'] - d) * vec[d - 1] for d in range(1, p['u'])) - (p['Q_real'] - n * p['u']) for p in per_u]
            assert min(slacks) >= 0 and sum(vec) == n * (n - 1) // 2
            signatures.append(dict(L=L, n=n, density=n/L, cap=max(vec), pair_sum=sum(vec),
                zeros=sum(x == 0 for x in vec), saturated=sum(x == max(vec) for x in vec),
                binding_windows=[p['u'] for p, s in zip(per_u, slacks) if s == 0],
                adjacent_equal=sum(a == b for a, b in zip(vec, vec[1:])),
                odd_distance_mass=sum(vec[::2]), even_distance_mass=sum(vec[1::2]), r=vec))
    assert len(rows) == 7 and len(signatures) == 5
    family = []
    for k in list(range(3, 41)) + [50, 64, 80, 100]:
        L, n, u = 4*k*k, 2*k*k, 2*k
        b = bounds(L, n)
        q = min_boundary_square(L, n, u)[0]
        assert 3*q == 12*k**4 - 4*k**3 + 3*k*k + k
        assert 6*(q - n*u) == (6*k*k - 5*k - 1)*u*(u-1)
        lo = k*k - (5*k+1)//6
        upper = k*k-k+2
        assert b['B_int'] <= upper and b['B_bd'] >= lo
        if k >= 14:
            assert b['B_bd'] > b['B_int']
        family.append(dict(k=k, **b, proved_integer_upper=upper, proved_boundary_lower=lo,
                           actual_gain=b['B_bd']-b['B_int'], proved_gain=-(-(k-1)//6)-2))
    for name, data in [('realizability_audit', audit), ('seven_case_hierarchy', rows),
                       ('seven_case_profiles', profiles), ('residual_signatures', signatures),
                       ('infinite_family_checks', family)]:
        (HERE / 'data' / (name + '.json')).write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    report = dict(status='PASS', seven_complete_cases=rows, independent_walk_windows=sum(len(r['per_u']) for r in profiles),
                  residual_signatures=signatures, infinite_family_numeric_checks=len(family),
                  seconds=round(time.perf_counter()-started, 3))
    (HERE/'tests/revision_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    run()
