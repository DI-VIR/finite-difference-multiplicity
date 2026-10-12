"""Cross-check the explicit-candidate construction on finite graph fixtures."""
from _bootstrap import report
from np_reduction import verify

if __name__=='__main__':
    count=verify()
    assert count==406
    report(graph_fixtures=count,scope='Sanity checks; the universal proof is in main Section 7')
