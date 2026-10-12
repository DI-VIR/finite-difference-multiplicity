"""Compatibility entry point for the preserved complexity implementation."""
from _source_loader import load_module,run_file

if __name__=='__main__':
    run_file('complexity/np_reduction.py')
else:
    load_module(__name__,'complexity/np_reduction.py')
