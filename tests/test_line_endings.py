"""Regression checks for Windows, POSIX and legacy line endings."""
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from regenerate_tables import canonical_text

with tempfile.TemporaryDirectory() as folder:
    paths = []
    for i, ending in enumerate((b'\n', b'\r\n', b'\r')):
        p = Path(folder) / str(i)
        p.write_bytes((b'\xef\xbb\xbf' if i == 1 else b'') + ending.join((b'alpha,beta', b'1,2')))
        paths.append(p)
    assert len({canonical_text(p) for p in paths}) == 1
    paths[0].write_bytes(b'alpha,beta\n1,3')
    assert canonical_text(paths[0]) != canonical_text(paths[1])
print('PASS: normalized line endings and BOM; content differences still fail.')
