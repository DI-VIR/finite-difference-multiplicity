"""Compatibility entry point for the preserved milp implementation."""
from _source_loader import load_module,run_file

if __name__=='__main__':
    run_file('milp/milp_difference_optimum.py')
else:
    load_module(__name__,'milp/milp_difference_optimum.py')
