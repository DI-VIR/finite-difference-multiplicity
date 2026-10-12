"""Independent Fraction calculations, coverage and upper-witness contracts."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import csv
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from integrated_baselines import BASELINES, compute, table_specs
from regenerate_tables import read_tables, table_cells


def read(name):
    return json.loads((ROOT / 'data' / name).read_text(encoding='utf-8-sig'))


def ceil(x):
    return -(-x.numerator // x.denominator)


def independent_classical(L, n):
    real, integer, xu = [], [], [F(0)]
    for u in range(2, L + 1):
        N, M = L + u - 1, n * u
        q, eta = divmod(M, N)
        real.append(F(M * M, N * u * (u - 1)) - F(n, u - 1))
        integer.append(F((N - eta) * q ** 2 + eta * (q + 1) ** 2 - M, u * (u - 1)))
        xu.extend((F(n * n, L + u) - F(n + 1, u), F(2 * (u * n - L), u * (u - 1)) - F(2, u)))
    return dict(B_ph=ceil(F(n * (n - 1), 2 * (L - 1))),
                B_CS=ceil(max(real)), B_int=ceil(max(integer)), B_Xu=ceil(max(xu)))


def main():
    rows, summary, protocol = compute()
    assert rows == read('integrated_baselines.json')
    assert summary == read('baseline_summary.json')
    assert protocol == read('baseline_protocol.json')
    old_xu = {(r['L'], r['n']): r for r in read('xu_comparison.json')}
    old_grid = {(r['L'], r['n']): r for r in read('bound_landscape.json')}
    old_literature = {(r['L'], r['n']): r for r in read('literature_grid_comparison.json')}
    exact = {(r['L'], r['n']): r for r in read('exact_map_nontrivial.json')}
    medium = {(r['L'], r['n']): r for r in read('verified_medium.json')}
    cells = read('published_diameters.json')
    assert len(cells) == 565
    assert len(rows) == len({(r['L'], r['n']) for r in rows}) == 702
    assert {(r['L'], r['n']) for r in rows} == set(old_grid) | set(exact) | set(medium) == set(old_xu)
    for row in rows:
        L, n = pair = (row['L'], row['n'])
        for key, value in independent_classical(L, n).items():
            assert row[key] == value, (pair, key)
        for key in ('B_Xu', 'B_CS', 'B_int', 'B_bd'):
            assert row[key] == old_xu[pair][key], (pair, key)
        table_candidates = [cell['g'] + 1 for cell in cells if cell['n'] <= n and cell['diameter'] >= L]
        table = max([0, *table_candidates])
        assert row['B_GO_table'] == table
        assert row['B_prior'] == max(table, 0, 2 * n - L - 1)
        assert row['B_base'] == max(row[key] for key in BASELINES[:-1])
        assert row['B_full'] == max(row['B_base'], row['B_bd'])
        assert max(row['B_bd'], table) == max(row['B_bd'], row['B_prior'])
        assert row['B_bd'] >= row['B_int'] >= row['B_CS']
        for cell in row['GO_driving_cells']:
            assert cell['g'] + 1 == table and cell['n'] <= n and cell['diameter'] >= L
        if pair in old_grid:
            for key in ('B_ph', 'B_CS', 'B_int', 'B_bd', 'u_bd'):
                assert row[key] == old_grid[pair][key]
            assert row['B_GO_table'] == old_literature[pair]['table_lower']
            assert row['B_GO_lemma'] == old_literature[pair]['GO_lemma']
        if pair in exact:
            assert row['B_base'] == row['verified_upper'] == exact[pair]['Delta']
        if pair in medium:
            witness = list(map(int, medium[pair]['witness'].split()))
            assert len(witness) == len(set(witness)) == n and min(witness) >= 0 and max(witness) < L
            differences = Counter(abs(a - b) for i, a in enumerate(witness) for b in witness[i + 1:])
            assert max(differences.values()) == row['verified_upper'] == medium[pair]['best_feasible']
    grid = summary['grid']['comparison']['B_base']
    assert (grid['stronger'], grid['tied'], grid['weaker']) == (184, 269, 9)
    assert summary['grid']['gain_histogram'] == {'-1': 9, '0': 269, '1': 180, '2': 4}
    assert summary['exact']['baseline_tight'] == 252
    assert summary['principal']['gain_histogram'] == {'0': 6, '1': 3}
    for domain in ('grid', 'exact', 'principal', 'unique_union'):
        subset = rows if domain == 'unique_union' else [r for r in rows if domain in r['domains']]
        for name in BASELINES:
            c = summary[domain]['comparison'][name]
            gains = [r['B_bd'] - r[name] for r in subset]
            assert c['stronger'] == sum(g > 0 for g in gains)
            assert c['tied'] == sum(g == 0 for g in gains)
            assert c['weaker'] == sum(g < 0 for g in gains)
            assert c['minimum_difference'] == min(gains) and c['maximum_difference'] == max(gains)
            assert c['stronger'] + c['tied'] + c['weaker'] == len(subset)
    tables = read_tables()
    expected = tuple(table_cells(spec) for spec in table_specs(rows, summary))
    assert (tables['8'], tables['S2b']) == expected
    with (ROOT / 'data/integrated_baselines.csv').open(encoding='utf-8-sig', newline='') as handle:
        exported = list(csv.DictReader(handle))
    assert len(exported) == len(rows)
    for saved, row in zip(exported, rows):
        for key, value in row.items():
            if isinstance(value, list):
                assert json.loads(saved[key]) == value
            elif value is None:
                assert saved[key] == ''
            else:
                assert saved[key] == str(value)
    print('PASS: 702 exact baseline comparisons, 565 source cells, nine upper witnesses, separate-domain counts and CSV/table contracts.')


if __name__ == '__main__':
    main()
