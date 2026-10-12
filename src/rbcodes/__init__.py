from __future__ import absolute_import
__version__ = "3.0.0"


def selftest(verbose=False):
    """Check that rbcodes is installed and working; returns True if all checks pass.
    Same as the ``rb_selftest`` command."""
    from rbcodes._selftest import selftest as _selftest
    return _selftest(verbose=verbose)
