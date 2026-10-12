"""Load preserved implementations while keeping historical flat commands valid."""
from pathlib import Path
import importlib.util
import runpy
import sys

ROOT=Path(__file__).resolve().parent

def load_module(name,relative):
    path=ROOT/relative
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module

def run_file(relative):
    runpy.run_path(str(ROOT/relative),run_name='__main__')
