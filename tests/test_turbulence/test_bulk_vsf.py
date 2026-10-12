"""Known-answer tests for rbcodes.turbulence.bulk and rbcodes.turbulence.vsf."""
import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from rbcodes.ifu.vmap import VelocityMap
from rbcodes.turbulence.bulk import (ArctanDisk, Constant, Plane, Quadratic, TwoLobe, arctan_disk,
                                     compare_bulk, two_lobe)
from rbcodes.turbulence.vsf import StructureFunction, _smoothed_powerlaw


def _amp(n, alpha, smooth=0.):
    k = np.hypot(*np.meshgrid(np.fft.fftfreq(n), np.fft.fftfreq(n)))
    amp = np.zeros_like(k)
    amp[k > 0] = k[k > 0] ** (-(alpha + 2) / 2)
    return amp * np.exp(-0.5 * (2 * np.pi * k * smooth) ** 2)     # Gaussian smoothing, sigma = smooth px


def grf(n, alpha, seed, rms=20., smooth=0.):
    """2-D Gaussian random field with P(k) ~ k^-(alpha+2) (periodic box), so S2 ~ r^alpha well
    inside the box and away from the pixel scale; optionally smoothed by a Gaussian."""
    rng = np.random.default_rng(seed)
    f = np.real(np.fft.ifft2(_amp(n, alpha, smooth) * np.fft.fft2(rng.normal(size=(n, n)))))
    return f / f.std() * rms


def expected_s2(sf, n, alpha, rms=20., smooth=0.):
    """Exact ensemble S2 per bin for the pairs of ``sf`` (positions in pixels), from the
    lattice power spectrum of ``grf`` (normalised to the unsmoothed variance rms^2)."""
    C = np.real(np.fft.ifft2(_amp(n, alpha, smooth) ** 2))
    C0 = np.real(np.fft.ifft2(_amp(n, alpha) ** 2))[0, 0]
    C = C / C0 * rms ** 2
    dx = np.rint(sf.x[sf.j] - sf.x[sf.i]).astype(int) % n
    dy = np.rint(sf.y[sf.j] - sf.y[sf.i]).astype(int) % n
    val = 2 * (C[0, 0] - C[dy, dx])
    g = sf.b >= 0
    num = np.bincount(sf.b[g], weights=val[g], minlength=sf.nb)
    cnt = np.bincount(sf.b[g], minlength=sf.nb)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(cnt > 0, num / cnt, np.nan)


def fit_exact(sf, S, rmin, rmax, **kw):
    """fit_powerlaw on an exact S2 (1% diagonal errors)."""
    S0, C0 = sf.S, sf.C
    sf.S, sf.C = S, np.diag((0.01 * np.nan_to_num(S, nan=1.)) ** 2 + 1e-12)
    try:
        return sf.fit_powerlaw(rmin, rmax, **kw)
    finally:
        sf.S, sf.C = S0, C0


def vmap(v, err=1., mask=None, kpp=None):
    v = np.asarray(v, float)
    e = np.broadcast_to(np.asarray(err, float), v.shape).copy()
    return VelocityMap(v, e, np.ones(v.shape, bool) if mask is None else mask, kpc_per_pix=kpp)


# ---------------------------------------------------------------------------
# bulk
# ---------------------------------------------------------------------------

class TestBulk:

    def test_polynomials_exact(self):
        yy, xx = np.mgrid[0:15, 0:20].astype(float)
        vm = vmap(30 + 4 * xx - 2 * yy + 0.1 * xx * yy, kpp=1.5)
        q = Quadratic()
        res = vm.subtract_bulk(q)
        assert np.max(np.abs(res.v)) < 1e-9
        p = vm.subtract_bulk(Plane())
        g, pa = p.meta['bulk']['p'][1:], None
        c = vm.subtract_bulk(Constant())
        assert abs(c.v[c.mask].mean()) < 1e-9
        assert Plane().fit(xx.ravel(), yy.ravel(), (30 + 4 * xx - 2 * yy).ravel()).gradient[0] == pytest.approx(np.hypot(4, 2))

    @pytest.mark.parametrize('kpp', [None, 1.49])
    def test_arctan_and_twolobe_recover_and_are_unit_invariant(self, kpp):
        yy, xx = np.mgrid[0:30, 0:36].astype(float)
        s = 1. if kpp is None else kpp
        true_d = [17 * s, 14 * s, 0.6, 0.6, 120., 4. * s, 10.]
        true_t = [18 * s, 15 * s, 0.4, 5., 80., 3. * s]
        for model, func, true in ((ArctanDisk(), arctan_disk, true_d), (TwoLobe(), two_lobe, true_t)):
            vm = vmap(func(true, xx * s, yy * s), kpp=kpp)
            res = vm.subtract_bulk(model)
            assert np.sqrt(np.mean(res.v ** 2)) < 0.5, type(model).__name__
            np.testing.assert_allclose(model.predict(xx * s, yy * s), func(true, xx * s, yy * s), atol=0.5)

    def test_compare_bulk_picks_the_right_model(self):
        rng = np.random.default_rng(0)
        yy, xx = np.mgrid[0:30, 0:36].astype(float)
        disk = vmap(arctan_disk([17, 14, 0.6, 0.6, 120., 4., 10.], xx, yy) + rng.normal(0, 3, xx.shape), err=3.)
        res, best = compare_bulk(disk)
        assert best == 'arctan' and res['plane']['dbic'] > 50
        plane = vmap(5 + 2 * xx - yy + rng.normal(0, 3, xx.shape), err=3.)
        res, best = compare_bulk(plane, dict(none=Constant(), plane=Plane(), quad=Quadratic()))
        assert best == 'plane' and res['plane']['rms'] == pytest.approx(3., rel=0.1)


# ---------------------------------------------------------------------------
# structure function
# ---------------------------------------------------------------------------

class TestStructureFunction:

    def test_matches_brute_force(self):
        rng = np.random.default_rng(1)
        v = rng.normal(0, 10, (9, 11))
        err = rng.uniform(1, 3, v.shape)
        mask = rng.random(v.shape) > 0.2
        vm = vmap(v, err, mask, kpp=2.)
        sf = StructureFunction(vm, n_bins=6, min_pairs=5)
        ys, xs = np.where(vm.mask)
        acc = {b: [] for b in range(6)}
        for a in range(len(ys)):
            for b in range(a + 1, len(ys)):
                r = 2. * np.hypot(xs[a] - xs[b], ys[a] - ys[b])
                k = np.searchsorted(sf.edges, r, side='right') - 1
                if 0 <= k < 6 and sf.ok[k]:
                    acc[k].append((v[ys[a], xs[a]] - v[ys[b], xs[b]]) ** 2 - err[ys[a], xs[a]] ** 2 - err[ys[b], xs[b]] ** 2)
        for k in range(6):
            if sf.ok[k]:
                assert sf.S[k] == pytest.approx(np.mean(acc[k]))

    def test_debiasing_removes_noise(self):
        rng = np.random.default_rng(2)
        err = np.linspace(5, 15, 30)[None, :] * np.ones((30, 1))
        vm = vmap(rng.normal(0, 1, err.shape) * err, err)
        sf = StructureFunction(vm, n_bins=10)
        g = sf.ok
        assert np.all(np.abs(sf.S[g]) < 3 * sf.err[g] + 5)
        np.testing.assert_allclose(sf.S_raw[g], 2 * np.mean(err ** 2), rtol=0.15)
        ivw = sf.s2('ivw')
        assert np.all(np.abs(ivw[g]) < 3 * sf.err[g] + 5)

    def test_s2_matches_expectation_for_random_fields(self):
        Ss, alphas = [], []
        for seed in range(8):
            rng = np.random.default_rng(50 + seed)
            vm = vmap(grf(96, 2 / 3, seed)[:48, :48] + rng.normal(0, 3, (48, 48)), err=3., kpp=1.)
            sf = StructureFunction(vm, n_bins=14)
            Ss.append(sf.S)
            alphas.append(sf.fit_powerlaw(2., 10.)['alpha'])
        Ss = np.array(Ss)
        exp = expected_s2(sf, 96, 2 / 3)
        g = sf.ok & (sf.rc < 20)
        mean, sem = Ss.mean(0)[g], Ss.std(0, ddof=1)[g] / np.sqrt(len(Ss))
        assert np.all(np.abs(mean - exp[g]) < 3 * sem + 0.02 * exp[g])
        assert np.mean(alphas) == pytest.approx(fit_exact(sf, exp, 2., 10.)['alpha'], abs=0.08)

    def test_jackknife_errors_have_the_right_size(self):
        # valid when the field's correlation length is much smaller than the region
        # (jackknife blocks cannot see variance from modes larger than the region)
        S, E = [], []
        for seed in range(30):
            rng = np.random.default_rng(100 + seed)
            f = gaussian_filter(rng.normal(0, 1, (40, 40)), 1.0, mode='wrap')
            vm = vmap(f / 0.2821 * 10 + rng.normal(0, 2, f.shape), err=2.)     # fixed scale (std = 1/(2 sqrt(pi) sigma))
            sf = StructureFunction(vm, n_bins=8)
            k = np.where(sf.ok)[0][3]
            S.append(sf.S[k])
            E.append(sf.err[k])
        assert 0.7 < np.std(S) / np.median(E) < 1.4

    def test_cross_estimator_has_no_noise_bias(self):
        rng = np.random.default_rng(3)
        field = grf(64, 1.0, 7)[:32, :32]
        a = vmap(field + rng.normal(0, 8, field.shape), 8.)
        b = vmap(field + rng.normal(0, 8, field.shape), 8.)
        truth = StructureFunction(vmap(field, 1e-6), n_bins=8, estimator='raw')
        cross = StructureFunction(a, n_bins=8, estimator='cross', other=b)
        g = cross.ok & truth.ok
        assert np.all(np.abs(cross.S[g] - truth.S[g]) < 3 * cross.err[g])
        assert np.all(cross.S_raw[g] > truth.S[g] + 64)
        with pytest.raises(ValueError):
            StructureFunction(a, estimator='cross')

    def test_directional(self):
        yy, xx = np.mgrid[0:30, 0:30].astype(float)
        sf = StructureFunction(vmap(3. * xx, 1e-3), n_bins=8)     # pair value 9 dx^2
        d = sf.directional(0., tol=np.pi / 12)
        Sp, Sq = d['parallel'][0], d['perpendicular'][0]
        g = np.isfinite(Sp) & np.isfinite(Sq)
        # parallel: dx^2 >= cos^2(15 deg) r^2, perpendicular: dx^2 <= sin^2(15 deg) r^2
        assert np.all(Sq[g] < 9 * np.sin(np.pi / 12) ** 2 * sf.edges[1:][g] ** 2)
        assert np.all(Sp[g] > 9 * np.cos(np.pi / 12) ** 2 * sf.edges[:-1][g] ** 2 * 0.99)
        d = sf.directional(np.pi / 2, tol=np.pi / 12)                # axis along y: roles swap
        np.testing.assert_allclose(d['parallel'][0][g], Sq[g])

    def test_source_plane_and_images(self):
        v = grf(64, 1.0, 9)[:20, :20]
        vm = vmap(v, 1., kpp=1.)
        img = np.where(vm.x < 10, 1, 2)
        vm.set_source_positions(0.5 * vm.x, 0.5 * vm.y, img)
        sfi = StructureFunction(vm, n_bins=8)
        sfs = StructureFunction(vm, n_bins=8, positions='source')
        assert sfs.pixel == pytest.approx(0.5) and sfs.r_geo == pytest.approx(0.5 * sfi.r_geo, rel=0.3)
        n1 = (img == 1).sum()
        assert len(sfs.i) == 2 * n1 * (n1 - 1) // 2                   # no pairs across images
        both = StructureFunction(vm, n_bins=8, positions='source', cross_images=True)
        assert len(sfi.i) == len(sfs.i) and len(both.i) == 400 * 399 // 2
        one = StructureFunction(vm.select(image=1), n_bins=8, positions='source')
        assert len(one.i) == n1 * (n1 - 1) // 2

    def test_noise_floor_fit(self):
        sf = StructureFunction(vmap(np.zeros((40, 40)), 1.), n_bins=12)
        S = 50. * (sf.rc / 10.) ** 0.8 + 72.
        f = fit_exact(sf, S, 1., 30., noise_floor=True)
        assert f['N'] == pytest.approx(72., rel=1e-3) and f['alpha'] == pytest.approx(0.8, abs=1e-3)
        assert fit_exact(sf, S, 1., 30.)['alpha'] < 0.5            # ignoring the floor flattens the slope

    def test_psf_model(self):
        r = np.array([0.5, 1., 3., 10., 30.])
        np.testing.assert_allclose(_smoothed_powerlaw(r, 0.8, 0.) / r ** 0.8,
                                   _smoothed_powerlaw(np.array([1.]), 0.8, 0.)[0], rtol=1e-6)
        # exact binned S2 of a power-law field smoothed by a Gaussian (sigma = 4 px): the PSF fit
        # recovers the input slope (to the accuracy of the lattice field, a few hundredths);
        # the plain fit is far off
        sf = StructureFunction(vmap(np.zeros((64, 64)), 1.), n_bins=16)
        for alpha in (0.67, 1.0):
            sm = expected_s2(sf, 512, alpha, smooth=4.0)
            plain = fit_exact(sf, sm, 6., 24.)['alpha']
            psf = fit_exact(sf, sm, 6., 24., psf_fwhm=4.0 * 2.3548)['alpha']
            assert abs(psf - alpha) < 0.06 and abs(plain - alpha) > 0.3, (alpha, psf, plain)

    def test_r_turn_and_increments(self):
        vm = vmap(grf(64, 1.0, 11)[:20, :20], 1.)
        sf = StructureFunction(vm, n_bins=10, min_pairs=1)
        S = np.minimum(sf.rc, 5.)
        assert sf.r_turn(S=S, frac=0.9) == pytest.approx(4.5, rel=0.05)
        assert len(sf.increments(1.)) == 2 * 20 * 19 + 2 * 19 * 19        # +-0.5 px: also diagonals (1.41)
