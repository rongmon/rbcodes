"""Tests for rbcodes.halo against analytic NFW results."""
import astropy.units as u
import numpy as np
import pytest
from astropy import constants as cnst

from rbcodes.halo.halo_profile import NFW_escape_vel, NFW_property_from_Mvir
from rbcodes.halo.rb_nfw import rb_nfw


def _m_nfw(r, rs, norm):
    x = r / rs
    return norm * (np.log(1. + x) - x / (1. + x))


@pytest.fixture(scope='module')
def nfw():
    return rb_nfw(m200=1e12, c=10., z=0.5)


class TestRbNfw:

    def test_enclosed_mass_matches_analytic(self, nfw):
        # M_r is integrated from the first grid point (1 kpc), not from 0
        r, r200 = nfw['r'], nfw['r200']
        rs = r200 / nfw['c']
        norm = nfw['m200'] / (np.log(1. + nfw['c']) - nfw['c'] / (1. + nfw['c']))
        expected = _m_nfw(r, rs, norm) - _m_nfw(r[0], rs, norm)
        np.testing.assert_allclose(nfw['M_r'][10:], expected[10:], rtol=5e-3)

    def test_mass_within_r200_is_m200(self, nfw):
        m_r200 = np.interp(nfw['r200'], nfw['r'], nfw['M_r'])
        assert m_r200 == pytest.approx(1e12, rel=0.01)

    def test_vcirc_peak_at_2p163_rs(self, nfw):
        rs = nfw['r200'] / nfw['c']
        dr = nfw['r'][1] - nfw['r'][0]
        assert abs(nfw['maxvcircr'][0] - 2.163 * rs) <= dr

    def test_outputs_finite(self, nfw):
        for key in ('rho', 'vcirc', 'M_r', 'sig_v', 'vesc'):
            assert np.all(np.isfinite(nfw[key])), key


class TestNFWEscapeVel:

    Mvir = 1e12 * u.solMass
    Rvir = 250. * u.kpc
    cvir = 10. * u.one

    def test_untruncated_matches_analytic(self):
        r = np.array([10., 50., 200.]) * u.kpc
        a = self.Rvir / self.cvir
        fofc = np.log(1 + 10.) - 10. / 11.
        phi = -cnst.G * self.Mvir / fofc * np.log(1 + r / a) / r
        expected = np.sqrt(2 * np.abs(phi)).to(u.km / u.s)
        got = NFW_escape_vel(r, self.Mvir, self.Rvir, self.cvir)
        np.testing.assert_allclose(got.value, expected.value, rtol=1e-10)

    def test_truncated_is_point_mass_outside_and_continuous(self):
        rtrunc = 2. * self.Rvir
        r = np.array([0.999, 1.001, 2., 3.]) * rtrunc
        v = NFW_escape_vel(r, self.Mvir, self.Rvir, self.cvir, truncated=2.)
        assert v[0].value == pytest.approx(v[1].value, rel=5e-3)
        # outside the truncation radius v_esc^2 r is constant (point mass)
        vr = (v[1:] ** 2 * r[1:]).to(u.km ** 2 / u.s ** 2 * u.kpc).value
        np.testing.assert_allclose(vr, vr[0], rtol=1e-10)

    def test_property_from_mvir_runs(self):
        out = NFW_property_from_Mvir(r=50. * u.kpc, Mvir=self.Mvir, z=0.5)
        assert out['rvir'].unit.is_equivalent(u.kpc)
        assert 50. < out['vesc'].to(u.km / u.s).value < 2000.
