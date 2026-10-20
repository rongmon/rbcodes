"""
Installation self-test for rbcodes.

Run from a terminal with ``rb_selftest`` or from Python with ``rbcodes.selftest()``.
It needs no display and no test dependencies, and takes a few seconds:

1. reports the versions of Python, rbcodes and the main dependencies;
2. imports every rbcodes module;
3. checks that every source file compiles without SyntaxWarning;
4. checks that every NumPy / SciPy / astropy / pandas / matplotlib attribute used
   in the source exists in the installed versions;
5. runs a few numerical checks against known answers.

The checks in 2-4 are also used by ``tests/test_package_health.py``.
"""
import ast
import importlib
import inspect
import os
import pkgutil
import platform
import sys
import tempfile
import time
import traceback
import warnings
from importlib import metadata
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
DEPENDENCIES = ('numpy', 'scipy', 'matplotlib', 'astropy', 'pandas', 'PyQt5', 'linetools',
                'emcee', 'photutils', 'corner', 'scikit-learn', 'tqdm', 'regions', 'ipywidgets')
CHECKED_PACKAGES = ('numpy', 'scipy', 'astropy', 'pandas', 'matplotlib')


def _headless():
    """GUI modules must import without a display."""
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    os.environ.setdefault('MPLBACKEND', 'Agg')


# ---------------------------------------------------------------------------
# Package-wide checks
# ---------------------------------------------------------------------------

def module_names():
    import rbcodes
    return sorted(m.name for m in pkgutil.walk_packages(rbcodes.__path__, 'rbcodes.')
                  if '.examples' not in m.name)


def source_files():
    return sorted(p for p in PACKAGE_DIR.rglob('*.py') if '__pycache__' not in p.parts)


def compile_without_syntax_warning(path):
    with warnings.catch_warnings():
        warnings.simplefilter('error', SyntaxWarning)
        compile(Path(path).read_text(encoding='utf-8'), str(path), 'exec')


def _checked(name):
    return name.split('.')[0] in CHECKED_PACKAGES


def _resolve(full):
    """Import ``full`` as a module, or as ``module.attribute``; None if neither exists."""
    try:
        return importlib.import_module(full)
    except ImportError:
        pass
    mod, _, attr = full.rpartition('.')
    try:
        return getattr(importlib.import_module(mod), attr)
    except (ImportError, AttributeError, ValueError):
        return None


def _aliases(tree):
    """Map local names to objects of the checked packages, plus unresolvable from-imports."""
    aliases, missing = {}, []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if _checked(a.name):
                    # 'import a.b' binds 'a'; 'import a.b as c' binds 'c' to a.b
                    target = a.name if a.asname else a.name.split('.')[0]
                    aliases[a.asname or target] = _resolve(target)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            if _checked(node.module):
                for a in node.names:
                    if a.name == '*':
                        continue
                    obj = _resolve(f'{node.module}.{a.name}')
                    if obj is None:
                        missing.append(f'{node.lineno}: from {node.module} import {a.name}')
                    aliases[a.asname or a.name] = obj
    return {k: v for k, v in aliases.items() if v is not None}, missing


def _dotted(node):
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
        return parts[::-1]
    return None


def _has(obj, attr):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if hasattr(obj, attr):
            return True
    if inspect.ismodule(obj):
        # submodules that are not imported by their parent package
        try:
            importlib.import_module(f'{obj.__name__}.{attr}')
            return True
        except ImportError:
            return False
    return False


def _get(obj, attr):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            return getattr(obj, attr)
        except AttributeError:
            return importlib.import_module(f'{obj.__name__}.{attr}')


def _owner(obj):
    return getattr(obj, '__name__' if inspect.ismodule(obj) else '__module__', '') or ''


def find_removed_attributes(path):
    """Attributes of numpy/scipy/astropy/pandas/matplotlib used in ``path`` that do not exist.

    Follows dotted names through modules and classes (``pd.DataFrame.append``);
    attributes of instances (``df.append``) cannot be resolved statically.
    """
    tree = ast.parse(Path(path).read_text(encoding='utf-8'))
    aliases, missing = _aliases(tree)
    missing = set(missing)
    if not aliases:
        return sorted(missing)
    # local names rebound by assignment are no longer the imported object
    rebound = {t.id for n in ast.walk(tree) if isinstance(n, (ast.Assign, ast.AugAssign, ast.AnnAssign))
               for t in (n.targets if isinstance(n, ast.Assign) else [n.target]) if isinstance(t, ast.Name)}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Attribute):
            continue
        parts = _dotted(node)
        if not parts or parts[0] not in aliases or parts[0] in rebound:
            continue
        obj = aliases[parts[0]]
        for attr in parts[1:]:
            if not (inspect.ismodule(obj) or inspect.isclass(obj)):
                break
            if not _has(obj, attr):
                if _checked(_owner(obj)):
                    missing.add(f'{node.lineno}: {".".join(parts)}')
                break
            obj = _get(obj, attr)
    return sorted(missing)


# ---------------------------------------------------------------------------
# Numerical checks against known answers
# ---------------------------------------------------------------------------

# Data files read at run time; they must be listed in setup.cfg [options.package_data]
# or they are missing from non-editable installs.
PACKAGE_DATA = ('halo/smhm_med_params.txt', 'IGM/lines/atom_full.dat', 'IGM/lines/lls.lst',
                'GUIs/multispecviewer/line_options.conf', 'GUIs/zgui/gui_linelists.ascii',
                'GUIs/interactive_continuum_fit_help/help.html',
                'GUIs/zfind/templates/marz', 'GUIs/zfind/templates/pca')


def check_package_data():
    missing = [f for f in PACKAGE_DATA if not (PACKAGE_DIR / f).exists()]
    assert not missing, f'missing data files: {missing}'
    for d in PACKAGE_DATA:
        if (PACKAGE_DIR / d).is_dir():
            assert list((PACKAGE_DIR / d).glob('*.fits')), f'no templates in {d}'


def check_line_list():
    from rbcodes.IGM.rb_setline import rb_setline
    out = rb_setline(2796.3, 'closest')
    assert out['name'] == 'MgII 2796' and abs(out['wave'] - 2796.352) < 1e-3, out


def check_equivalent_width():
    import numpy as np
    from rbcodes.IGM.compute_EW import compute_EW
    wrest = 1215.67
    lam = np.linspace(1213., 1218.5, 2000)
    v = 2.998e5 * (lam - wrest) / wrest
    flux = 1. - 0.5 * np.exp(-(v / 30.) ** 2)
    expected = np.trapezoid(1. - flux, x=lam)
    out = compute_EW(lam, flux, wrest, [-200., 200.], 0.01 * np.ones_like(lam), f0=0.4164)
    assert abs(out['ew_tot'] / expected - 1.) < 0.01, (out['ew_tot'], expected)


def check_nfw_mass():
    import numpy as np
    from rbcodes.halo.rb_nfw import rb_nfw
    nfw = rb_nfw(m200=1e12, c=10., z=0.5)
    m_r200 = np.interp(nfw['r200'], nfw['r'], nfw['M_r'])
    assert abs(m_r200 / 1e12 - 1.) < 0.01, m_r200


def angular_diameter_distance_12(cosmo, z1, z2):
    """astropy D_A between two redshifts (the method name changed in astropy 8)."""
    try:
        return cosmo.angular_diameter_distance(z1, z2)
    except TypeError:
        return cosmo.angular_diameter_distance_z1z2(z1, z2)


def check_cosmology():
    from astropy.cosmology import FlatLambdaCDM
    from rbcodes.lensing.lens_ang_sep import ang_D12, cosmic_D
    cosmo = FlatLambdaCDM(H0=70., Om0=0.3, Tcmb0=0.)
    mpc = 3.0856775814913673e22
    da = cosmic_D(0.3, 0.7, 1.5)['D_A'] / mpc
    d12 = ang_D12(0.3, 0.7, 0.5, 1.5) / mpc
    assert abs(da / cosmo.angular_diameter_distance(1.5).value - 1.) < 2e-3, da
    assert abs(d12 / angular_diameter_distance_12(cosmo, 0.5, 1.5).value - 1.) < 2e-3, d12


def check_lens_raytrace():
    """Singular isothermal sphere: a point on the Einstein ring maps to the lens centre."""
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    from rbcodes.lensing.lens_ang_sep import model_delens
    n, pix, theta_e = 101, 0.1, 2.0
    c = (n - 1) / 2
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [150., 2.]
    w.wcs.crpix = [c + 1, c + 1]
    w.wcs.cdelt = [-pix / 3600., pix / 3600.]
    y, x = np.mgrid[0:n, 0:n] - c
    r = np.hypot(x, y)
    r[r == 0] = 1.
    with tempfile.TemporaryDirectory() as d:
        fx, fy = os.path.join(d, 'dplx.fits'), os.path.join(d, 'dply.fits')
        fits.PrimaryHDU(theta_e * x / r, header=w.to_header()).writeto(fx)
        fits.PrimaryHDU(theta_e * y / r, header=w.to_header()).writeto(fy)
        s = model_delens(fx, fy, zl=0.5, zs_o=2.0)
        ra, dec = w.pixel_to_world_values(np.array([c + theta_e / pix]), np.array([c]))
        sra, sdec = s.raytrace_new_z(ra, dec, zs_n=2.0)
        bx, by = w.world_to_pixel_values(sra, sdec)
    assert abs(bx[0] - c) < 1e-6 and abs(by[0] - c) < 1e-6, (bx, by)


def check_ifu_cube_io():
    import numpy as np
    from astropy.io import fits
    from rbcodes.ifu.io import load_fits
    nw, ny, nx = 50, 6, 7
    h = fits.Header()
    h['INSTRUME'] = 'MUSE'
    for k, v in dict(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CTYPE3='AWAV', CUNIT3='Angstrom',
                     CRVAL1=150., CRVAL2=2., CRVAL3=4750., CRPIX1=1., CRPIX2=1., CRPIX3=1.,
                     CDELT1=-0.2 / 3600, CDELT2=0.2 / 3600, CDELT3=1.25).items():
        h[k] = v
    flux = np.random.default_rng(0).normal(size=(nw, ny, nx))
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 'cube.fits')
        fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(flux, header=h, name='DATA'),
                      fits.ImageHDU(np.ones_like(flux), header=h, name='STAT')]).writeto(f)
        cube = load_fits(f)
    assert cube.flux.shape == (nw, ny, nx) and cube.var is not None
    assert abs(cube.wave[0] - 4750.) < 1e-9 and abs(cube.wave[-1] - (4750. + 1.25 * (nw - 1))) < 1e-9
    assert cube.wcs is not None and cube.wcs.naxis == 2


def check_ifu_line_fit():
    import numpy as np
    from rbcodes.ifu.lines import LineModel, lsf_muse
    from rbcodes.ifu.measure import GaussianFitter
    lm = LineModel('OII', z=0.785, windows=dict(fit=(6620., 6690.)), air=False, lsf=lsf_muse)
    wave = np.arange(6620., 6690., 1.25)
    true = (0.02, 4., 1.2, -75., 35.)          # continuum within the default +-5 x noise bound
    r = GaussianFitter(lm).fit_spectrum(wave, lm(wave, *true), np.full(wave.size, 0.01))
    assert r['status'] == 0 and np.allclose(r['p'], true, rtol=1e-4, atol=1e-4), r.get('p')


NUMERICAL_CHECKS = (check_package_data, check_line_list, check_equivalent_width, check_nfw_mass,
                    check_cosmology, check_lens_raytrace, check_ifu_cube_io, check_ifu_line_fit)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def _versions():
    import rbcodes
    lines = [f'Python {platform.python_version()} ({sys.executable})',
             f'rbcodes {rbcodes.__version__} ({PACKAGE_DIR})']
    for dep in DEPENDENCIES:
        try:
            lines.append(f'{dep} {metadata.version(dep)}')
        except metadata.PackageNotFoundError:
            lines.append(f'{dep} NOT INSTALLED')
    return lines


def selftest(verbose=False):
    """
    Check that rbcodes is installed and working. Prints a report and returns
    True if every check passed.

    verbose: also print the full traceback of each failure.
    """
    _headless()
    t0 = time.time()
    print('rbcodes self-test')
    print('-' * 60)
    for line in _versions():
        print('  ' + line)
    print('-' * 60)

    failures = []

    def record(name, problems):
        status = 'PASS' if not problems else 'FAIL'
        print(f'{status}  {name}')
        for p in problems:
            print(f'        {p}')
        if problems:
            failures.append(name)

    problems = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for name in module_names():
            try:
                importlib.import_module(name)
            except Exception as e:
                problems.append(f'{name}: {type(e).__name__}: {e}')
    record('import all modules', problems)

    problems = []
    for path in source_files():
        try:
            compile_without_syntax_warning(path)
        except SyntaxError as e:
            problems.append(f'{path.relative_to(PACKAGE_DIR)}: {e}')
    record('source compiles without warnings', problems)

    problems = []
    for path in source_files():
        try:
            problems += [f'{path.relative_to(PACKAGE_DIR)}:{m}' for m in find_removed_attributes(path)]
        except SyntaxError:
            pass  # reported above
    record('no removed numpy/scipy/astropy/pandas/matplotlib APIs', problems)

    for check in NUMERICAL_CHECKS:
        name = check.__name__[len('check_'):].replace('_', ' ')
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                check()
            record(name, [])
        except Exception as e:
            msg = [f'{type(e).__name__}: {e}']
            if verbose:
                msg += traceback.format_exc().rstrip().split('\n')
            record(name, msg)

    print('-' * 60)
    if failures:
        print(f'FAILED: {len(failures)} check(s) failed ({time.time() - t0:.1f} s). '
              'Please include this full report when reporting the problem.')
    else:
        print(f'ALL CHECKS PASSED ({time.time() - t0:.1f} s)')
    return not failures


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Check that rbcodes is installed and working.')
    parser.add_argument('-v', '--verbose', action='store_true', help='show full tracebacks')
    args = parser.parse_args()
    sys.exit(0 if selftest(verbose=args.verbose) else 1)


if __name__ == '__main__':
    main()
