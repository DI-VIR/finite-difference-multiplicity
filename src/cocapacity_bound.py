"""Compatibility entry point for the preserved bounds implementation."""
from _source_loader import load_module,run_file

if __name__=='__main__':
    run_file('bounds/cocapacity_bound.py')
else:
    load_module(__name__,'bounds/cocapacity_bound.py')
