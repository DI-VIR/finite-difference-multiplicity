"""Re-export frozen CSV tables without a separate catalog or LaTeX files.

This is a deterministic export, not a substitute for underlying numerical checks.
"""
from pathlib import Path
import csv

HERE=Path(__file__).resolve().parents[1]


def canonical_text(path):
    """Ignore a UTF-8 BOM and platform line endings, but preserve all content."""
    return Path(path).read_text(encoding='utf-8-sig').replace('\r\n', '\n').replace('\r', '\n')


def table_cells(spec):
    """Use CSV's text representation for comparisons with computed table rows."""
    def cell(value):
        return '' if value is None else str(value)
    return dict(headers=[cell(value) for value in spec['headers']],
                rows=[[cell(value) for value in row] for row in spec['rows']])


def read_tables(folder=None):
    folder = Path(folder) if folder else HERE/'tables'
    tables = {}
    for path in sorted(folder.glob('Table_*.csv')):
        with path.open(encoding='utf-8-sig', newline='') as handle:
            rows = list(csv.reader(handle))
        assert rows and all(len(row) == len(rows[0]) for row in rows), path.name
        tables[path.stem.removeprefix('Table_')] = dict(headers=rows[0], rows=rows[1:])
    assert tables, 'No frozen CSV tables found'
    return tables


def write_tables(tables, output=None):
    out = Path(output) if output else HERE/'tables'
    out.mkdir(parents=True, exist_ok=True)
    for key, spec in tables.items():
        cells = table_cells(spec)
        assert all(len(row) == len(cells['headers']) for row in cells['rows']), key
        with (out/f'Table_{key}.csv').open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.writer(handle, lineterminator='\n')
            writer.writerow(cells['headers'])
            writer.writerows(cells['rows'])
    return len(tables)


def regenerate(output=None, verify=True):
    out = Path(output) if output else HERE/'recomputed/tables'
    tables = read_tables()
    write_tables(tables, out)
    if verify:
        for key in tables:
            name = f'Table_{key}.csv'
            assert canonical_text(out/name) == canonical_text(HERE/'tables'/name), key
    return len(tables)


if __name__=='__main__':print('Canonical table exports:',regenerate())
