"""Compatibility entry point for the preserved bounds implementation."""
from _source_loader import load_module,run_file

if __name__=='__main__':
    run_file('bounds/joint_verification.py')
else:
    load_module(__name__,'bounds/joint_verification.py')
