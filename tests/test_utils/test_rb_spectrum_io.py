"""Reading the example 1D spectra with rb_spectrum."""
import contextlib
import gc
import io
import warnings
from pathlib import Path

import numpy as np
import pytest

from rbcodes.utils.rb_spectrum import rb_spectrum

DATA = Path(__file__).resolve().parents[2] / 'example-data'
SPECTRA = ['test.fits', 'test_norm.fits', 'test_normalized.fits', 'galaxy1.fits', 'sdss1.fits',
           'sdss2.fits', 'spec-0398-51789-0282.fits', 'J075103+424211a_F.fits',
           'jw01243001001_01101_00001_nrca5_x1d.fits',
           'eiger-mock-data/stacked_1D_14099_visit1_mod_a.fits']


@pytest.mark.parametrize('name', SPECTRA)
def test_reads_example_spectrum_and_closes_file(name):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', ResourceWarning)
        with contextlib.redirect_stdout(io.StringIO()):
            sp = rb_spectrum.from_file(str(DATA / name))
        gc.collect()
    assert not getattr(sp, '_read_failed', False)
    wave, flux = np.asarray(sp.wavelength), np.asarray(sp.flux)
    assert wave.size > 10 and wave.shape == flux.shape
    # some files pad the wavelength array with NaN (as stored in the file)
    good = np.isfinite(wave)
    assert good.sum() > 0.5 * wave.size and np.all(np.diff(wave[good]) != 0)
    assert np.isfinite(flux).sum() > 0.5 * flux.size
    leaked = [w for w in caught if issubclass(w.category, ResourceWarning) and name in str(w.message)]
    assert not leaked, leaked[0].message
