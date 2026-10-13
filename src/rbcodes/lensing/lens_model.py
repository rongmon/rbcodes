"""
Lens models for velocity maps: source-plane positions, magnification and multiple images.

A thin layer over ``lens_ang_sep.model_delens`` (which does the ray tracing). It adds
lens-model realizations, ray tracing of a ``VelocityMap`` into source-plane coordinates
(kpc), the magnification tensor, and matching of spaxels across multiple images.

    from rbcodes.lensing.lens_model import LensModel, match_images
    lens = LensModel('dplx.fits', 'dply.fits', zl=0.43, zs=2.76)            # best fit only
    lens = LensModel.from_directory('Lens_Model/range', zl=0.43, zs=2.76)  # best fit + realizations
    vm_src = lens.raytrace(vm, z=1.70, images={1: 'arc1.reg', 2: 'arc2.reg'})
    for k, vm_k in lens.raytrace_realizations(vm, z=1.70):                  # k = 0 is the best fit
        ...
    pairs = match_images(vm_src, tol=0.5)                                   # same-position test

Ray tracing moves positions only: velocities and errors stay with their image-plane spaxels.
Astrometry of the cube and the lens model must already agree (it is not adjusted here).

Source-plane coordinates: offsets from a fixed sky position ``centre`` (default: the mean
source-plane position of the kept spaxels under the best-fit model), x = -(RA - RA0) cos(Dec0), y = Dec - Dec0
(x increases westward, like the pixel x axis of a north-up, east-left image), converted to
proper kpc at the plane redshift with a flat-LambdaCDM cosmology (H0 = 70, the model's Omega_m),
matching the distances used by the ray tracer.
"""
import glob
import os
import re

import numpy as np

from .lens_ang_sep import model_delens

__all__ = ['LensModel', 'match_images']


class LensModel:
    """
    Parameters
    ----------
    fits_x, fits_y : str
        Best-fit deflection maps (arcsec, along the pixel axes, with a celestial WCS).
    zl : float
        Lens redshift.
    zs : float
        Source redshift the deflection maps were computed for.
    realizations : list of (fits_x, fits_y), optional
        Other lens-model realizations (e.g. MCMC samples). Without them only the best fit is used.
    Omega_m, Omega_l : float
        Cosmology of the lens model (default 0.3, 0.7).
    """

    def __init__(self, fits_x, fits_y, zl, zs, realizations=None, Omega_m=0.3, Omega_l=0.7):
        self.zl, self.zs = float(zl), float(zs)
        self.Omega_m, self.Omega_l = Omega_m, Omega_l
        self.files = [(fits_x, fits_y)] + [tuple(r) for r in (realizations or [])]
        self._models = {}

    @classmethod
    def from_fits(cls, fits_x, fits_y, zl, zs, **kw):
        return cls(fits_x, fits_y, zl, zs, **kw)

    @classmethod
    def from_directory(cls, path, zl, zs, pattern_x='dplx_*.fits', pattern_y='dply_*.fits', best='1', **kw):
        """
        Best fit and realizations from a folder of numbered maps (e.g. ``dplx_1.fits`` ...
        ``dplx_121.fits``); ``best`` is the label of the best-fit model.
        """
        def label(f, pat):
            rx = re.escape(pat).replace(r'\*', '(.+)')
            m = re.fullmatch(rx, os.path.basename(f))
            return m.group(1) if m else None
        xs = {label(f, pattern_x): f for f in glob.glob(os.path.join(path, pattern_x))}
        ys = {label(f, pattern_y): f for f in glob.glob(os.path.join(path, pattern_y))}
        labels = sorted(set(xs) & set(ys), key=lambda s: (len(s), s))
        if best not in labels:
            raise FileNotFoundError(f"No best-fit maps '{pattern_x.replace('*', best)}' / "
                                    f"'{pattern_y.replace('*', best)}' in {path}.")
        others = [(xs[k], ys[k]) for k in labels if k != best]
        return cls(xs[best], ys[best], zl, zs, realizations=others, **kw)

    @property
    def n_realizations(self):
        """Number of lens models (best fit + realizations)."""
        return len(self.files)

    def model(self, realization=0):
        """The ``model_delens`` object of a realization (0 = best fit)."""
        if realization not in self._models:
            fx, fy = self.files[realization]
            m = model_delens(fx, fy, zl=self.zl, zs_o=self.zs)
            m.update_cosmology(Omega_m=self.Omega_m, Omega_l=self.Omega_l)
            self._models[realization] = m
        return self._models[realization]

    # --- positions --------------------------------------------------------------

    def source_positions(self, ra, dec, z, realization=0, interp='nearest'):
        """Source-plane RA, Dec (deg) at redshift z of image-plane positions (NaN outside the maps)."""
        return self.model(realization).raytrace_new_z(ra, dec, zs_n=z, interp=interp)

    def kpc_per_arcsec(self, z):
        from astropy.cosmology import LambdaCDM
        cosmo = LambdaCDM(H0=70., Om0=self.Omega_m, Ode0=self.Omega_l, Tcmb0=0.)
        return cosmo.kpc_proper_per_arcmin(z).value / 60.

    def raytrace(self, vm, z, realization=0, images=None, centre=None, interp='nearest'):
        """
        New VelocityMap with source-plane positions (kpc) at redshift ``z`` for every spaxel.

        images : {id: bool (ny, nx) mask or ds9 region file/text}, optional
            Which multiple image each spaxel belongs to (others get -1).
        centre : (ra, dec) in degrees, optional
            Origin of the source-plane coordinates (default ``source_centre``: the mean source-plane
            position of the kept spaxels under the best-fit model, the same for every realization).
        """
        if vm.wcs is None:
            raise ValueError('The VelocityMap has no WCS; it is needed to ray-trace its spaxels.')
        ra, dec = vm.wcs.pixel_to_world_values(vm.x, vm.y)
        sra, sdec = self.source_positions(ra.ravel(), dec.ravel(), z, realization, interp)
        sra, sdec = np.reshape(sra, vm.shape), np.reshape(sdec, vm.shape)
        if centre is None:
            centre = self.source_centre(vm, z, interp) if realization != 0 else _mean_pos(sra, sdec, vm.mask)
        k = self.kpc_per_arcsec(z) * 3600.
        bx = -(sra - centre[0]) * np.cos(np.radians(centre[1])) * k
        by = (sdec - centre[1]) * k
        image = _image_ids(vm, images) if images is not None else vm.image
        new = vm.copy()
        new.set_source_positions(bx, by, image)
        new.meta['lens'] = dict(zl=self.zl, zs=self.zs, z_plane=float(z), realization=int(realization),
                                files=list(self.files[realization]), centre=list(centre), interp=interp)
        new.provenance.append(f'ray traced to z={z} with lens model {realization} '
                              f'({os.path.basename(self.files[realization][0])})')
        return new

    def raytrace_realizations(self, vm, z, images=None, centre=None, interp='nearest'):
        """Yield (k, VelocityMap) for every lens model, k = 0 the best fit; same centre for all."""
        if centre is None:
            centre = self.source_centre(vm, z, interp)
        for k in range(self.n_realizations):
            yield k, self.raytrace(vm, z, k, images, centre, interp)

    def source_centre(self, vm, z, interp='nearest'):
        """Mean source-plane (RA, Dec) of the kept spaxels of ``vm`` under the best-fit model (realization 0)."""
        ra, dec = vm.wcs.pixel_to_world_values(vm.x[vm.mask], vm.y[vm.mask])
        sra, sdec = self.source_positions(ra, dec, z, 0, interp)
        return _mean_pos(sra, sdec, np.ones(np.shape(sra), bool))

    # --- magnification ----------------------------------------------------------

    def magnification(self, ra, dec, z, realization=0, interp='nearest'):
        """
        Lensing Jacobian A = I - s dalpha/dtheta and magnification mu = 1/det(A) at image-plane
        positions, for a source at redshift z (s = deflection_scale(z)); derivatives by central
        differences of the deflection maps, read at the nearest map pixel (``interp='nearest'``)
        or interpolated (``'bilinear'``; better near critical curves, where mu changes fast).

        Returns dict(a11, a12, a21, a22, mu, kappa, gamma) — kappa = 1 - (a11 + a22)/2,
        gamma = sqrt(((a11 - a22)/2)^2 + ((a12 + a21)/2)^2).
        """
        m = self.model(realization)
        px, py = m.wcs_x.world_to_pixel_values(np.asarray(ra, float), np.asarray(dec, float))
        px, py = np.asarray(px), np.asarray(py)
        ny, nx = m.x_def.shape
        if interp == 'nearest':
            ix = np.rint(np.nan_to_num(px, nan=-1)).astype(int)
            iy = np.rint(np.nan_to_num(py, nan=-1)).astype(int)
            ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
            look = lambda g: g[iy[ok], ix[ok]]
        elif interp == 'bilinear':
            from scipy.ndimage import map_coordinates
            ok = np.isfinite(px) & np.isfinite(py) & (px >= 0) & (px <= nx - 1) & (py >= 0) & (py <= ny - 1)
            look = lambda g: map_coordinates(g, np.vstack([py[ok], px[ok]]), order=1, mode='nearest')
        else:
            raise ValueError(f"interp must be 'nearest' or 'bilinear', got {interp!r}")
        if not hasattr(m, '_grad'):
            sx, sy = m.pix_scale_x * 3600., m.pix_scale_y * 3600.       # arcsec per pixel
            daxdy, daxdx = np.gradient(m.x_def)
            daydy, daydx = np.gradient(m.y_def)
            m._grad = (daxdx / sx, daxdy / sy, daydx / sx, daydy / sy)
        s = m.deflection_scale(z)
        out = {}
        for name, g, one in (('a11', m._grad[0], 1.), ('a12', m._grad[1], 0.),
                             ('a21', m._grad[2], 0.), ('a22', m._grad[3], 1.)):
            a = np.full(px.shape, np.nan)
            a[ok] = one - s * look(g)
            out[name] = a
        det = out['a11'] * out['a22'] - out['a12'] * out['a21']
        with np.errstate(divide='ignore', invalid='ignore'):
            out['mu'] = 1. / det
        out['kappa'] = 1. - 0.5 * (out['a11'] + out['a22'])
        out['gamma'] = np.hypot(0.5 * (out['a11'] - out['a22']), 0.5 * (out['a12'] + out['a21']))
        return out


def _mean_pos(ra, dec, mask):
    ok = mask & np.isfinite(ra) & np.isfinite(dec)
    if not ok.any():
        raise ValueError('No kept spaxel falls inside the deflection maps.')
    return float(np.mean(ra[ok])), float(np.mean(dec[ok]))


def _image_ids(vm, images):
    ids = np.full(vm.shape, -1, int)
    for k, reg in images.items():
        if isinstance(reg, np.ndarray):
            m = reg.astype(bool)
        else:
            from rbcodes.ifu.noise import sky_mask_from_region
            m = sky_mask_from_region(reg, vm.wcs, vm.shape)
        if (ids[m] >= 0).any():
            raise ValueError(f'Image {k} overlaps another image.')
        ids[m] = int(k)
    return ids


def match_images(vm, tol):
    """
    Pairs of kept spaxels in different images whose source-plane positions are within ``tol`` kpc.

    Returns dict with arrays: ``a``, ``b`` (flat indices into the map), ``image_a``, ``image_b``,
    ``dist`` (kpc), ``dv`` = v_a - v_b, ``err`` = sqrt(sigma_a^2 + sigma_b^2), ``pull`` = dv/err.
    Spaxels that map to the same source position should agree within ``err``.
    """
    from scipy.spatial import cKDTree
    if not vm.has_source_positions or vm.image is None:
        raise ValueError('Needs source-plane positions and image ids (LensModel.raytrace with images=).')
    keep = vm.mask & np.isfinite(vm.beta_x) & np.isfinite(vm.beta_y) & (vm.image >= 0)
    idx = np.flatnonzero(keep)
    pts = np.c_[vm.beta_x.ravel()[idx], vm.beta_y.ravel()[idx]]
    img = vm.image.ravel()[idx]
    pairs = np.array(sorted(cKDTree(pts).query_pairs(tol)), dtype=int).reshape(-1, 2)
    pairs = pairs[img[pairs[:, 0]] != img[pairs[:, 1]]]
    a, b = idx[pairs[:, 0]], idx[pairs[:, 1]]
    v, e = vm.v.ravel(), vm.v_err.ravel()
    dv, err = v[a] - v[b], np.hypot(e[a], e[b])
    return dict(a=a, b=b, image_a=vm.image.ravel()[a], image_b=vm.image.ravel()[b],
                dist=np.hypot(*(pts[pairs[:, 0]] - pts[pairs[:, 1]]).T), dv=dv, err=err, pull=dv / err)
