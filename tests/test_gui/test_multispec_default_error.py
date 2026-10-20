"""rb_multispec: spectra without an error array get the default error, from files and from data."""
import numpy as np
import pytest
from astropy.io import fits

pytest.importorskip('PyQt5')
from PyQt5.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope='module')
def app():
    return QApplication.instance() or QApplication([])


def _arrays():
    w = np.linspace(4000., 5000., 300)
    return w, 1. + 0.1 * np.sin(w / 50.)


def test_from_data_gets_default_error(app):
    from rbcodes.utils.rb_spectrum import rb_spectrum
    from rbcodes.GUIs.multispecviewer.rb_multispec import from_data
    from rbcodes.GUIs.multispecviewer.io_manager import DEFAULT_ERROR_FRACTION
    w, f = _arrays()
    no_err = rb_spectrum(wavelength=w, flux=f, filename='ifu_extraction', meta={'ra': 51.86})
    with_err = rb_spectrum(wavelength=w, flux=f, error=np.full_like(f, 0.2), filename='has_error')
    with pytest.warns(UserWarning, match='no error array'):
        win = from_data([no_err, with_err], show=False)
    a, b = win.spectra
    assert DEFAULT_ERROR_FRACTION == 0.05
    np.testing.assert_allclose(a.sig.value, 0.05 * f)
    assert a.meta.get('ra') == 51.86 and a.filename == 'ifu_extraction'      # metadata kept
    np.testing.assert_allclose(b.sig.value, 0.2)                               # existing errors untouched


def test_file_gets_default_error(app, tmp_path):
    from rbcodes.GUIs.multispecviewer.io_manager import IOManager
    w, f = _arrays()
    h = fits.Header()
    h['CRVAL1'], h['CDELT1'], h['CRPIX1'] = w[0], w[1] - w[0], 1.
    p = tmp_path / 'noerr.fits'
    fits.PrimaryHDU(f, header=h).writeto(p)
    spectra, err = IOManager().load_fits_files([str(p)])
    assert spectra and spectra[0].sig is not None
    np.testing.assert_allclose(spectra[0].sig.value, 0.05 * spectra[0].flux.value)
