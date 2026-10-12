"""Check that every distributed module imports without the master workspace."""
from pathlib import Path
import importlib
import json
import sys

folder=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
sys.path.insert(0,str(folder/'src'))
modules=[]
for path in sorted((folder/'src').glob('*.py')):
    importlib.import_module(path.stem)
    modules.append(path.stem)
result=dict(status='PASS',modules=modules,package=folder.name)
canonical=[]
for path in sorted((folder/'src').rglob('*.py')):
    relative=path.relative_to(folder/'src')
    if len(relative.parts)>1:
        name='.'.join(relative.with_suffix('').parts)
        module=importlib.import_module(name)
        assert Path(module.__file__).resolve()==path.resolve()
        canonical.append(relative.as_posix())
result['canonical_domain_modules']=canonical
if (folder/'archived_five_input/toy_state_recovery.py').exists():
    sys.path.insert(0,str(folder/'archived_five_input'))
    import toy_state_recovery
    toy_state_recovery.selftest()
    result['archived_five_input_selftest']='PASS'
print(json.dumps(result,indent=2))

