"""Check the risk-closure identities and weighted binary formulation."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from risk_closure_audit import audit

print(json.dumps(audit(), indent=2))
