"""Final submission details are independently reproducible from the scientific archive."""
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from final_details_audit import audit

print(json.dumps(audit(),indent=2))
