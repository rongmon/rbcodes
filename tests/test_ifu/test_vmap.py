"""Tests for rbcodes.ifu.vmap.VelocityMap."""
import numpy as np
import pytest
from astropy.wcs import WCS

from rbcodes.ifu.lines import LineModel, lsf_muse
from rbcodes.ifu.measure import GaussianFitter, MomentMeasurer, QualityMask
from rbcodes.ifu.vmap import VelocityMap

NY, NX = 10, 12


def _wcs():
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [231.94, 6.87]
    w.wcs.crpix = [6., 5.]
    w.wcs.cdelt = [-0.2 / 3600, 0.2 / 3600]
    return w


def _vm(seed=0):
    rng = np.random.default_rng(seed)
    v = rng.normal(0, 20, (NY, NX))
    mask = np.ones((NY, NX), bool)
    mask[0, 0] = False
    return VelocityMap(v, np.full((NY, NX), 5.), mask, width=np.full((NY, NX), 40.), kpc_per_pix=1.5,
                       psf_fwhm_kpc=5.2, line='OII', z=0.785, wcs=_wcs(), provenance=['test'])


class Plane:
    """Minimal bulk model with the fit/predict protocol."""
    def fit(self, x, y, v, ve):
        A = np.c_[np.ones_like(x), x, y]
        self.p = np.linalg.lstsq(A / ve[:, None], v / ve, rcond=None)[0]

    def predict(self, x, y):
        return self.p[0] + self.p[1] * x + self.p[2] * y

    def describe(self):
        return dict(kind='plane', p=self.p.tolist())


def test_from_fit_and_moments():
    lm = LineModel('OII', 0.785273, dict(fit=(6640., 6680.)), lsf=lsf_muse)
    wave = np.arange(6600., 6720., 1.25)
    rng = np.random.default_rng(1)
    flux = np.zeros((len(wave), 6, 7))
    dom = lm.fit_mask(wave)
    for iy in range(6):
        for ix in range(7):
            flux[dom, iy, ix] = lm(wave[dom], 0., 3., 1.2, 10. * ix, 40.)
    flux[:, :, 0] = 0.
    flux += rng.normal(0, 0.05, flux.shape)
    sig = np.full_like(flux, 0.05)
    fit = GaussianFitter(lm).fit_cube(flux, sig, wave=wave, verbose=False)
    q = QualityMask(min_island=5).apply(fit)
    vm = VelocityMap.from_fit(fit, q, kpc_per_pix=1.49)
    assert vm.mask.sum() == q.mask.sum() == 36
    np.testing.assert_array_equal(vm.code, q.code)
    np.testing.assert_array_equal(vm.v, fit['v'])                 # values kept, nothing filled
    assert vm.line == 'OII' and vm.z == pytest.approx(0.785273) and 'Gaussian fit' in vm.provenance[0]
    mom = MomentMeasurer(lm, vwin=(-250, 250)).measure_cube(flux, sig, wave=wave)
    vmm = VelocityMap.from_moments(mom, QualityMask(min_island=5).apply(mom))
    assert vmm.meta['method'] == 'moments' and vmm.mask.sum() > 30


def test_points_and_source_positions():
    vm = _vm()
    x, y, v, ve = vm.points()
    assert len(x) == NY * NX - 1 and x.max() == pytest.approx((NX - 1) * 1.5)
    with pytest.raises(ValueError):
        vm.points('source')
    bx, by = vm.x * 0.3, vm.y * 0.3
    bx[1, 1] = np.nan
    img = np.where(vm.x < 6, 1, 2)
    vm.set_source_positions(bx, by, img)
    xs, ys, vs, _ = vm.points('source')
    assert len(xs) == NY * NX - 2
    one = vm.select(image=1)
    assert one.mask.sum() == (img == 1).sum() - 1 and vm.mask.sum() == NY * NX - 1   # (0,0) is masked


def test_select_region_mask_and_ds9():
    vm = _vm()
    box = np.zeros((NY, NX), bool)
    box[2:5, 3:8] = True
    a = vm.select(box, name='box')
    b = vm.select('image\nbox(6,4,5,3,0)')                    # 1-based centre (6,4) -> rows 2-4, cols 3-7
    np.testing.assert_array_equal(a.mask, box)
    np.testing.assert_array_equal(b.mask, box)
    assert 'region box' in a.provenance[-1]


def test_subtract_bulk_recovers_plane():
    vm = _vm()
    vm.v = vm.v * 0 + 30. + 4. * vm.x * 1.5 - 2. * vm.y * 1.5 + np.random.default_rng(2).normal(0, 1, vm.shape)
    res = vm.subtract_bulk(Plane())
    np.testing.assert_allclose(res.meta['bulk']['p'], [30., 4., -2.], atol=0.3)
    assert abs(np.nanmean(res.v[res.mask])) < 0.3 and 'Plane' in res.provenance[-1]
    assert np.allclose(vm.v[vm.mask].mean(), 30. + 4. * vm.x[vm.mask].mean() * 1.5 - 2. * vm.y[vm.mask].mean() * 1.5, atol=0.5)


def test_save_load_roundtrip(tmp_path):
    vm = _vm()
    vm.set_source_positions(vm.x * 0.3, vm.y * 0.3, np.ones(vm.shape, int))
    vm = vm.select(image=1)
    p = tmp_path / 'vm.fits'
    vm.save(p)
    back = VelocityMap.load(p)
    for k in ('v', 'v_err', 'width', 'width_err', 'flux', 'flux_err', 'code', 'mask', 'beta_x', 'beta_y', 'image'):
        np.testing.assert_array_equal(getattr(back, k), getattr(vm, k), err_msg=k)
    assert (back.kpc_per_pix, back.psf_fwhm_kpc, back.line, back.z) == (1.5, 5.2, 'OII', 0.785)
    assert back.provenance == vm.provenance
    np.testing.assert_allclose(back.wcs.wcs.crval, [231.94, 6.87])


def test_plot_runs():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    _vm().plot(ax)
    _vm().plot(ax, 'code', show_masked=True, colorbar=False)
    plt.close(fig)
