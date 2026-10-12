"""
Package-wide health checks.

These catch breakage from new Python / library releases in code that has no
dedicated tests (GUIs, scripts, rarely used helpers):

* every module imports;
* every source file compiles without SyntaxWarning (e.g. invalid escapes);
* every attribute of NumPy, SciPy, astropy, pandas and matplotlib used in the
  source exists in the installed version (e.g. removed ``np.trapz``,
  ``pd.DataFrame.append``), and every ``from <pkg> import <name>`` resolves.

The checks live in ``rbcodes._selftest`` so that users can run them with
``rb_selftest`` after installing.
"""
import importlib
import textwrap

import pytest

from rbcodes import _selftest as st

PACKAGE_DIR = st.PACKAGE_DIR


def _rel(p):
    return str(p.relative_to(PACKAGE_DIR))


@pytest.mark.parametrize('name', st.module_names())
def test_module_imports(name):
    importlib.import_module(name)


@pytest.mark.parametrize('path', st.source_files(), ids=_rel)
def test_compiles_without_syntax_warning(path):
    st.compile_without_syntax_warning(path)


@pytest.mark.parametrize('path', st.source_files(), ids=_rel)
def test_library_attributes_exist(path):
    missing = st.find_removed_attributes(path)
    assert not missing, 'Removed or unknown library attributes:\n' + '\n'.join(missing)


def test_attribute_check_flags_removed_apis(tmp_path):
    """The checker must flag removed APIs and accept valid uses.

    Some APIs (np.trapz, cumtrapz, cwt) were only deprecated in the oldest supported
    versions, so the expected set depends on what is installed.
    """
    import numpy as np
    import scipy.integrate
    import scipy.signal
    src = tmp_path / 'probe.py'
    src.write_text(textwrap.dedent('''
        import numpy as np
        import pandas as pd
        import matplotlib as mpl
        from matplotlib import cm
        from astropy.io import fits
        from scipy import integrate
        from scipy.signal import cwt
        a = np.trapz([1, 2])
        b = pd.DataFrame.append
        c = cm.get_cmap('viridis')
        d = mpl.cm.get_cmap('viridis')
        e = fits.new_table([])
        f = integrate.cumtrapz([1, 2])
        g = np.not_a_numpy_function
        ok = (np.trapezoid, pd.concat, mpl.colormaps, fits.open, integrate.cumulative_trapezoid)
    '''))
    expected = {'pd.DataFrame.append', 'cm.get_cmap', 'mpl.cm.get_cmap', 'fits.new_table',
                'np.not_a_numpy_function'}
    if not hasattr(np, 'trapz'):
        expected.add('np.trapz')
    if not hasattr(scipy.integrate, 'cumtrapz'):
        expected.add('integrate.cumtrapz')
    if not hasattr(scipy.signal, 'cwt'):
        expected.add('from scipy.signal import cwt')
    flagged = {m.split(': ', 1)[1] for m in st.find_removed_attributes(src)}
    assert flagged == expected


@pytest.mark.parametrize('check', st.NUMERICAL_CHECKS, ids=lambda c: c.__name__)
def test_numerical_check(check):
    check()


def test_selftest_passes(capsys):
    assert st.selftest()
    assert 'ALL CHECKS PASSED' in capsys.readouterr().out
