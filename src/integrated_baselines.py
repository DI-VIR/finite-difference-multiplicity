"""Exact, same-target baseline comparison; never rerun historical searches."""
from pathlib import Path
from fractions import Fraction
import argparse
import csv
import json

from extremal import bounds

HERE = Path(__file__).resolve().parents[1]
BASELINES = ('B_ph', 'B_CS', 'B_int', 'B_Xu', 'B_prior', 'B_base')
GO_URL = 'https://arxiv.org/html/2605.14229v1'
XU_URL = 'https://arxiv.org/html/1706.05969v5'


def read(name):
    return json.loads((HERE / 'data' / name).read_text(encoding='utf-8-sig'))


def ceiling(value):
    return -(-value.numerator // value.denominator)


def xu_bound(length, cardinality):
    best, threshold, source = 0, 2, 'zero floor'
    for d in range(2, length + 1):
        candidates = (
            ('Theorem 3.1', Fraction(cardinality ** 2, length + d) - Fraction(cardinality + 1, d)),
            ('Proposition 2.2', Fraction(2 * (d * cardinality - length), d * (d - 1)) - Fraction(2, d)))
        for label, value in candidates:
            if ceiling(value) > best:
                best, threshold, source = ceiling(value), d, label
    return best, threshold, source


def prior_bound(length, cardinality, cells):
    eligible = [cell for cell in cells if cell['n'] <= cardinality]
    excluded = [cell for cell in eligible if cell['diameter'] > length - 1]
    table = max((cell['g'] + 1 for cell in excluded), default=0)
    driving = [cell for cell in excluded if cell['g'] + 1 == table]
    lemma = max(0, 2 * cardinality - length - 1)
    return dict(B_GO_lemma=lemma, B_GO_table=table, B_prior=max(lemma, table),
                GO_direct_cardinality_covered=any(cell['n'] == cardinality for cell in cells),
                GO_driving_cells=[dict(g=cell['g'], n=cell['n'], diameter=cell['diameter'])
                                  for cell in driving])


def summarize(rows):
    comparisons = {}
    for name in BASELINES:
        gains = [row['B_bd'] - row[name] for row in rows]
        comparisons[name] = dict(stronger=sum(g > 0 for g in gains),
                                 tied=sum(g == 0 for g in gains),
                                 weaker=sum(g < 0 for g in gains),
                                 minimum_difference=min(gains), maximum_difference=max(gains))
    return dict(pairs=len(rows), comparison=comparisons,
                reverse_pairs=[dict(L=r['L'], n=r['n'], B_bd=r['B_bd'], B_base=r['B_base'],
                                    baseline_winners=r['baseline_winners'])
                               for r in rows if r['B_bd'] < r['B_base']],
                gain_histogram={str(g): sum(r['B_bd'] - r['B_base'] == g for r in rows)
                                for g in sorted({r['B_bd'] - r['B_base'] for r in rows})},
                boundary_tight=sum(r['verified_upper'] == r['B_bd'] for r in rows),
                baseline_tight=sum(r['verified_upper'] == r['B_base'] for r in rows),
                combined_tight=sum(r['verified_upper'] == r['B_full'] for r in rows))


def compute():
    grid = read('bound_landscape.json')
    exact = read('exact_map_nontrivial.json')
    medium = read('verified_medium.json')
    cells = read('published_diameters.json')
    populations = {name: {(r['L'], r['n']) for r in rows}
                   for name, rows in (('grid', grid), ('exact', exact), ('principal', medium))}
    exact_index = {(r['L'], r['n']): r for r in exact}
    medium_index = {(r['L'], r['n']): r for r in medium}
    rows = []
    for length, cardinality in sorted(set.union(*populations.values())):
        row = bounds(length, cardinality)
        value, threshold, source = xu_bound(length, cardinality)
        row.update(B_Xu=value, D_Xu=threshold, Xu_source=source)
        row.update(prior_bound(length, cardinality, cells))
        row['B_base'] = max(row[name] for name in BASELINES[:-1])
        row['baseline_winners'] = [name for name in BASELINES[:-1] if row[name] == row['B_base']]
        row['B_full'] = max(row['B_base'], row['B_bd'])
        row['net_boundary_gain'] = row['B_bd'] - row['B_base']
        pair = (length, cardinality)
        row['domains'] = [name for name, pairs in populations.items() if pair in pairs]
        row['verified_upper'] = None
        row['upper_source'] = None
        if pair in medium_index:
            row['verified_upper'] = medium_index[pair]['best_feasible']
            row['upper_source'] = 'verified_medium.json (historical checked witness)'
        if pair in exact_index:
            optimum = exact_index[pair]['Delta']
            if row['verified_upper'] is not None:
                assert row['verified_upper'] == optimum
            row['verified_upper'] = optimum
            row['upper_source'] = 'exact_map_nontrivial.json (exhaustive anchored enumeration)'
        if row['verified_upper'] is not None:
            assert row['B_full'] <= row['verified_upper'], row
        rows.append(row)
    summary = {name: summarize([r for r in rows if name in r['domains']]) for name in populations}
    summary['unique_union'] = summarize(rows)
    summary['scope'] = ('Exact deterministic comparisons on three existing domains, not independent '
                        'random observations or new historical solver/search runs. Domains overlap.')
    protocol = dict(
        target='Delta(L,n), positive differences, exactly n marks in positions 0,...,L-1',
        arithmetic='Integer and Fraction arithmetic; ceilings precede maximization; all integer windows 2,...,L.',
        grid_source='bound_landscape.json (462 prespecified pairs)',
        exact_source='exact_map_nontrivial.json (252 nontrivial pairs)',
        principal_source='verified_medium.json (nine principal pairs)',
        unique_pairs=len(rows), published_GO_cells=len(cells),
        sources=dict(
            B_ph='ceil(n(n-1)/(2(L-1))); unordered-pair pigeonhole',
            B_CS='max_u ceil(((nu)^2/(L+u-1)-nu)/(u(u-1))); optimized classical window bound',
            B_int='max_u ceil((Nq^2+(2q+1)eta-nu)/(u(u-1))); N=L+u-1, nu=Nq+eta',
            B_bd='max_u ceil((Q_bd(L,n,u)-nu)/(u(u-1))); main Theorem 1',
            B_Xu=dict(reference=21, url=XU_URL, version='v5',
                      formulas=['max_D ceil(n^2/(L+D)-(n+1)/D)',
                                'max_D ceil(2(Dn-L)/(D(D-1))-2/D)'],
                      thresholds='2<=D<=L; maximum also includes zero',
                      theorem_numbers=['Theorem 3.1', 'Proposition 2.2']),
            B_prior=dict(reference=2, url=GO_URL, version='v1',
                         lemma='max(0,2n-L-1), Lemma 6',
                         table='max({0} union {g+1: stored G(g,k)>L-1, k<=n}), Table 1',
                         combination='max(GO lemma, GO finite-table inversion)',
                         no_exclusion='A zero table value is lack of a certificate, not Delta=0.'),
            B_base='max(B_ph,B_CS,B_int,B_Xu,B_prior); strongest of these declared baselines only',
            B_full='max(B_base,B_bd); full combined lower certificate',
            B_comb='Historical max(B_bd,B_GO_table) is preserved; no experiment stopping target is changed.'),
        restrictions=['No claimed domination of all results of a cited paper or of unevaluated LP/SDP models.',
                      'Known exact diameter exclusions are prior information, not newly discovered values.',
                      'GO direct-cardinality coverage and monotone-transfer-only strata remain distinct.',
                      'A verified upper witness is not an optimum unless it meets a valid lower bound.',
                      'No finite constant is assigned to an unspecified asymptotic remainder.'],
        primary_source_check_date='2026-10-10')
    return rows, summary, protocol


def table_specs(rows, summary):
    principal = [row for row in rows if 'principal' in row['domains']]
    order = [(r['L'], r['n']) for r in read('verified_medium.json')]
    indexed = {(r['L'], r['n']): r for r in principal}
    columns = ('B_ph', 'B_CS', 'B_int', 'B_Xu', 'B_prior', 'B_bd', 'verified_upper')
    main = dict(
        caption=('Same-target lower bounds and independently verified upper witnesses. B_prior is the '
                 'maximum of [2, Lemma 6] and the Table 1 inversion; B_Xu uses the two finite inequalities '
                 'of [21] specified in Section S17. Classical bounds are optimized over all integer '
                 'windows. An upper witness establishes an optimum only when it meets a valid lower bound.'),
        headers=['$(L,n)$', 'B_ph', 'B_CS', 'B_int', 'B_Xu', 'B_prior', 'B_bd', 'Verified upper'],
        rows=[[f'({p[0]},{p[1]})', *[indexed[p][c] for c in columns]] for p in order])
    labels = dict(B_ph='B_ph', B_CS='B_CS', B_int='B_int', B_Xu='B_Xu', B_prior='B_prior', B_base='B_base')
    si = dict(
        caption=('Boundary bound versus each declared comparator on the fixed 462-pair grid. '
                 'Stronger/tied/weaker use the exact integer difference B_bd minus the comparator. '
                 'B_base is the pointwise maximum of the five preceding baselines. The GO coverage '
                 'strata and all reverse source cells remain in Tables S2--S3.'),
        headers=['Comparator', 'Stronger', 'Tied', 'Weaker', 'Minimum difference', 'Maximum difference'],
        rows=[[labels[name], *[summary['grid']['comparison'][name][field] for field in
                               ('stronger', 'tied', 'weaker', 'minimum_difference', 'maximum_difference')]]
              for name in BASELINES])
    return main, si


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Export new comparison records and two table specifications.')
    args = parser.parse_args()
    rows, summary, protocol = compute()
    if args.write:
        write_json(HERE / 'data/integrated_baselines.json', rows)
        write_json(HERE / 'data/baseline_summary.json', summary)
        write_json(HERE / 'data/baseline_protocol.json', protocol)
        with (HERE / 'data/integrated_baselines.csv').open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
            writer.writeheader()
            for row in rows:
                writer.writerow({key: json.dumps(value) if isinstance(value, list) else value
                                 for key, value in row.items()})
        from regenerate_tables import write_tables
        main_table, comparison_table = table_specs(rows, summary)
        write_tables({'8': main_table, 'S2b': comparison_table})
    else:
        assert rows == read('integrated_baselines.json')
        assert summary == read('baseline_summary.json')
        assert protocol == read('baseline_protocol.json')
        from regenerate_tables import read_tables, table_cells
        tables = read_tables()
        expected = tuple(table_cells(spec) for spec in table_specs(rows, summary))
        assert (tables['8'], tables['S2b']) == expected
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
