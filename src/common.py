"""Package-local data paths. No dependency on the private master workspace."""
from pathlib import Path
import csv
import json

HERE = Path(__file__).resolve().parents[1]


def csv_read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        return list(csv.DictReader(handle))


def save(name, rows, protocol=None):
    out = HERE/'data'; out.mkdir(exist_ok=True)
    (out/(name+'.json')).write_text(json.dumps(rows,indent=2),encoding='utf-8')
    if rows:
        keys=list(dict.fromkeys(k for row in rows for k in row))
        with (out/(name+'.csv')).open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
            writer.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)
    if protocol:
        (out/(name+'_protocol.json')).write_text(json.dumps(protocol,indent=2),encoding='utf-8')
