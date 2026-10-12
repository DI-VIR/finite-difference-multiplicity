"""Compatibility entry point for the preserved exact implementation."""
from _source_loader import load_module,run_file

if __name__=='__main__':
    run_file('exact/extremal.py')
else:
    load_module(__name__,'exact/extremal.py')
