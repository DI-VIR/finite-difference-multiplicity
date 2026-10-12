"""Package-local test imports; no master-project path is used."""
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))

def load(name):
    return json.loads((ROOT/'data'/(name+'.json')).read_text(encoding='utf-8'))

def report(**items):
    print(json.dumps(dict(status='PASS',**items),indent=2))
