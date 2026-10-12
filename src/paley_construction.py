"""Compatibility entry point for the preserved constructions implementation."""
from _source_loader import load_module,run_file

if __name__=='__main__':
    run_file('constructions/paley_construction.py')
else:
    load_module(__name__,'constructions/paley_construction.py')
