"""rb_ifuview load options (wavelength range, air -> vacuum)."""
import numpy as np
import pytest
from astropy.io import fits

pytest.importorskip('PyQt5')
from PyQt5.QtWidgets import QApplication  # noqa: E402

from rbcodes.GUIs.ifuviewer import load_options  # noqa: E402


@pytest.fixture(scope='module')
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def cube_file(tmp_path):
    h = fits.Header()
    for k, v in dict(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CTYPE3='AWAV', CUNIT3='Angstrom', CRVAL1=150.,
                     CRVAL2=2., CRVAL3=5000., CRPIX1=1., CRPIX2=1., CRPIX3=1., CDELT1=-5e-5, CDELT2=5e-5,
                     CDELT3=1.25).items():
        h[k] = v
    p = tmp_path / 'c.fits'
    fits.HDUList([fits.PrimaryHDU(header=fits.Header({'INSTRUME': 'MUSE'})),
                  fits.ImageHDU(np.ones((400, 5, 6), np.float32), h, name='DATA'),
                  fits.ImageHDU(np.ones((400, 5, 6), np.float32), h, name='STAT')]).writeto(p)
    return str(p)


def test_small_files_are_not_asked(app, cube_file):
    assert load_options.ask_load_options(cube_file) == {}
    wmin, wmax, airvac, is_cube = load_options.cube_header_info(cube_file)
    assert (wmin, airvac, is_cube) == (5000., 'air', True) and wmax == pytest.approx(5000 + 399 * 1.25)


def test_dialog_defaults_and_options(app, cube_file, monkeypatch):
    monkeypatch.setattr(load_options, 'LARGE_FILE_BYTES', 10)        # treat the test file as large
    d = load_options.LoadOptionsDialog(cube_file)
    assert d.part.isChecked() and d.vac.isEnabled()
    d.lo.setValue(5100.)
    d.hi.setValue(5150.)
    d.vac.setChecked(True)
    assert d.options() == dict(wave_range=(5100., 5150.), air_to_vacuum=True)
    d.full.setChecked(True)
    assert d.options()['wave_range'] is None


def test_sidebar_loads_a_slab(app, cube_file):
    from rbcodes.GUIs.ifuviewer.sidebar import DatasetSidebar
    sb = DatasetSidebar()
    sb.add_dataset(cube_file, wave_range=(5100., 5150.), air_to_vacuum=True)
    c = sb.active_cube
    assert c.flux.shape == (41, 5, 6) and c.airvac == 'vac' and c.wave[0] > 5100.
