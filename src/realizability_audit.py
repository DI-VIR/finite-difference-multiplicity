"""Compatibility entry point for the preserved realizability implementation."""
from _source_loader import load_module,run_file

if __name__=='__main__':
    run_file('realizability/realizability_audit.py')
else:
    load_module(__name__,'realizability/realizability_audit.py')
