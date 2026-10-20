"""
VelocityMap: the velocity field of one line, with its errors, mask and positions.

Data and mocks, lensed or unlensed, all become a ``VelocityMap``; everything downstream
(bulk-flow removal, structure functions, summaries, inference) only uses ``VelocityMap``.

    from rbcodes.ifu.vmap import VelocityMap
    vm = VelocityMap.from_fit(fit, qm, kpc_per_pix=1.49, psf_fwhm_kpc=5.2)
    vm.save('OII_vmap.fits');  vm = VelocityMap.load('OII_vmap.fits')
    inner = vm.select(region='inner.reg')            # ds9 region, bool mask, or image=...
    x, y, v, v_err = vm.points()                     # kept spaxels, positions in kpc (or pixels)

Positions: image-plane pixel coordinates (0-based x = column, y = row) always; source-plane
coordinates beta (kpc) and an image id when a lens model has been applied
(``set_source_positions``). Rejected spaxels keep their values but are outside ``mask``;
nothing is filled or interpolated.
"""
import json

import numpy as np

__all__ = ['VelocityMap']

_ARRAYS = ('v', 'v_err', 'width', 'width_err', 'flux', 'flux_err', 'code')


class VelocityMap:
    """
    Parameters
    ----------
    v, v_err : 2-D arrays (km/s)
    mask : bool 2-D array, spaxels to use
    width, width_err, flux, flux_err, code : 2-D arrays, optional
    kpc_per_pix : float, optional
        Image-plane scale; without it image positions are in pixels.
    psf_fwhm_kpc : float, optional
    line, z : metadata
    wcs : astropy.wcs.WCS (2-D), optional
    provenance : list of str
    meta : dict of other settings (e.g. the fit settings)
    """

    def __init__(self, v, v_err, mask, width=None, width_err=None, flux=None, flux_err=None, code=None,
                 kpc_per_pix=None, psf_fwhm_kpc=None, line='', z=None, wcs=None, provenance=None, meta=None):
        v = np.asarray(v, dtype=float)
        shape = v.shape
        nan = np.full(shape, np.nan)
        self.v, self.v_err = v, np.asarray(v_err, dtype=float)
        self.width = nan.copy() if width is None else np.asarray(width, dtype=float)
        self.width_err = nan.copy() if width_err is None else np.asarray(width_err, dtype=float)
        self.flux = nan.copy() if flux is None else np.asarray(flux, dtype=float)
        self.flux_err = nan.copy() if flux_err is None else np.asarray(flux_err, dtype=float)
        self.code = np.full(shape, -1, int) if code is None else np.asarray(code).astype(int)
        self.mask = np.asarray(mask, dtype=bool) & np.isfinite(v) & np.isfinite(self.v_err)
        for name in _ARRAYS + ('mask',):
            if getattr(self, name).shape != shape:
                raise ValueError(f'{name} has shape {getattr(self, name).shape}, expected {shape}.')
        self.y, self.x = np.mgrid[0:shape[0], 0:shape[1]].astype(float)
        self.beta_x = self.beta_y = None
        self.image = None
        self.kpc_per_pix = kpc_per_pix
        self.psf_fwhm_kpc = psf_fwhm_kpc
        self.line, self.z, self.wcs = line, z, wcs
        self.provenance = list(provenance or [])
        self.meta = dict(meta or {})

    # --- construction -------------------------------------------------------

    @classmethod
    def _from_result(cls, res, quality, method, **kw):
        mask = quality.mask if hasattr(quality, 'mask') else np.asarray(quality, bool)
        code = getattr(quality, 'code', None)
        meta = dict(res.meta)
        if hasattr(quality, 'settings'):
            meta['quality'] = quality.settings
        prov = list(kw.pop('provenance', None) or [])
        prov.insert(0, f"{method} {meta.get('line', '')}: {int(mask.sum())} spaxels kept")
        return cls(res['v'], res['v_err'], mask, width=res['sig'], width_err=res['sig_err'],
                   flux=res['F'], flux_err=res['F_err'], code=code, line=meta.get('line', ''),
                   z=meta.get('z'), wcs=res.wcs, provenance=prov, meta=meta, **kw)

    @classmethod
    def from_fit(cls, fit, quality, kpc_per_pix=None, psf_fwhm_kpc=None, provenance=None):
        """From a ``GaussianFitter`` FitResult and a ``QualityMask`` result (or a bool mask).
        ``width`` is the intrinsic dispersion when the fit included the LSF."""
        return cls._from_result(fit, quality, 'Gaussian fit', kpc_per_pix=kpc_per_pix,
                                psf_fwhm_kpc=psf_fwhm_kpc, provenance=provenance)

    @classmethod
    def from_moments(cls, moments, quality, kpc_per_pix=None, psf_fwhm_kpc=None, provenance=None):
        """From a ``MomentMeasurer`` FitResult; ``width`` includes the instrumental width."""
        return cls._from_result(moments, quality, 'moments', kpc_per_pix=kpc_per_pix,
                                psf_fwhm_kpc=psf_fwhm_kpc, provenance=provenance)

    def copy(self):
        new = VelocityMap(self.v.copy(), self.v_err.copy(), self.mask.copy(), self.width.copy(),
                          self.width_err.copy(), self.flux.copy(), self.flux_err.copy(), self.code.copy(),
                          self.kpc_per_pix, self.psf_fwhm_kpc, self.line, self.z, self.wcs,
                          list(self.provenance), json.loads(json.dumps(self.meta, default=_jd)))
        if self.beta_x is not None:
            new.set_source_positions(self.beta_x.copy(), self.beta_y.copy(),
                                     None if self.image is None else self.image.copy())
        return new

    # --- positions ----------------------------------------------------------

    @property
    def shape(self):
        return self.v.shape

    @property
    def has_source_positions(self):
        return self.beta_x is not None

    def set_source_positions(self, beta_x, beta_y, image=None):
        """Source-plane coordinates (kpc, 2-D like v; NaN where unknown) and image ids (int, -1 = none)."""
        bx, by = np.asarray(beta_x, dtype=float), np.asarray(beta_y, dtype=float)
        if bx.shape != self.shape or by.shape != self.shape:
            raise ValueError('beta arrays must have the shape of the map.')
        self.beta_x, self.beta_y = bx, by
        self.image = None if image is None else np.asarray(image).astype(int)
        if self.image is not None and self.image.shape != self.shape:
            raise ValueError('image ids must have the shape of the map.')

    def points(self, positions='image'):
        """
        Kept spaxels as flat arrays (x, y, v, v_err).

        positions='image': pixel coordinates times kpc_per_pix (pixels if kpc_per_pix is None);
        positions='source': beta (kpc); spaxels without beta are dropped.
        """
        m = self.mask.copy()
        if positions == 'image':
            s = 1. if self.kpc_per_pix is None else self.kpc_per_pix
            x, y = self.x * s, self.y * s
        elif positions == 'source':
            if not self.has_source_positions:
                raise ValueError('No source-plane positions; apply a lens model first.')
            x, y = self.beta_x, self.beta_y
            m &= np.isfinite(x) & np.isfinite(y)
        else:
            raise ValueError("positions must be 'image' or 'source'.")
        return x[m], y[m], self.v[m], self.v_err[m]

    # --- selection and bulk flows -----------------------------------------

    def select(self, region=None, image=None, name=None):
        """
        New map keeping only spaxels inside ``region`` and/or of lensed image ``image``.

        region : bool (ny, nx) array, or a ds9 region file / text (needs ``wcs`` for sky regions).
        """
        new = self.copy()
        keep = np.ones(self.shape, bool)
        what = []
        if region is not None:
            if isinstance(region, np.ndarray):
                r = region.astype(bool)
                if r.shape != self.shape:
                    raise ValueError('region mask must have the shape of the map.')
            else:
                from .noise import sky_mask_from_region
                r = sky_mask_from_region(region, self.wcs, self.shape)
            keep &= r
            what.append(f'region {name or ("mask" if isinstance(region, np.ndarray) else region)}')
        if image is not None:
            if self.image is None:
                raise ValueError('No image ids; apply a lens model first.')
            keep &= self.image == image
            what.append(f'image {image}')
        new.mask &= keep
        new.provenance.append(f"select {', '.join(what)}: {int(new.mask.sum())} spaxels")
        return new

    def subtract_bulk(self, model, positions='image'):
        """
        New map with a bulk-flow model fitted to the kept spaxels and subtracted.

        ``model`` provides ``fit(x, y, v, v_err)`` and ``predict(x, y)`` (positions as in
        ``points``); e.g. the models in ``rbcodes.turbulence.bulk``. The fitted model is stored in
        ``meta['bulk']`` when it has ``describe()``.
        """
        x, y, v, ve = self.points(positions)
        model.fit(x, y, v, ve)
        if positions == 'image':
            s = 1. if self.kpc_per_pix is None else self.kpc_per_pix
            X, Y = self.x * s, self.y * s
        else:
            X, Y = self.beta_x, self.beta_y
        new = self.copy()
        new.v = self.v - model.predict(X, Y)
        new.mask &= np.isfinite(new.v)
        name = type(model).__name__
        if hasattr(model, 'describe'):
            new.meta['bulk'] = model.describe()
        new.provenance.append(f'subtract bulk {name} ({positions} plane)')
        return new

    # --- I/O ----------------------------------------------------------------

    def save(self, path, overwrite=True):
        """Multi-extension FITS: one image per array (2-D WCS on each), metadata as JSON."""
        from astropy.io import fits
        ph = fits.Header()
        ph['ORIGIN'] = 'rbcodes.ifu.vmap'
        ph['LINE'] = self.line
        if self.z is not None:
            ph['Z'] = self.z
        if self.kpc_per_pix is not None:
            ph['KPCPIX'] = (self.kpc_per_pix, 'kpc per pixel (image plane)')
        if self.psf_fwhm_kpc is not None:
            ph['PSFKPC'] = (self.psf_fwhm_kpc, 'PSF FWHM [kpc]')
        ph['PROVEN'] = json.dumps(self.provenance)
        ph['META'] = json.dumps(self.meta, default=_jd)
        wh = self.wcs.to_header() if self.wcs is not None else fits.Header()
        arrays = {k: getattr(self, k) for k in _ARRAYS}
        arrays['mask'] = self.mask.astype(np.int16)
        if self.has_source_positions:
            arrays['beta_x'], arrays['beta_y'] = self.beta_x, self.beta_y
            if self.image is not None:
                arrays['image'] = self.image
        hdus = [fits.PrimaryHDU(header=ph)]
        for k, a in arrays.items():
            h = wh.copy()
            h['EXTNAME'] = k.upper()
            hdus.append(fits.ImageHDU(a, header=h))
        fits.HDUList(hdus).writeto(path, overwrite=overwrite)

    @classmethod
    def load(cls, path):
        from astropy.io import fits
        from astropy.wcs import WCS
        with fits.open(path, memmap=False) as hdul:
            ph = hdul[0].header
            d = {h.name.lower(): h.data for h in hdul[1:]}
            wcs = WCS(hdul[1].header, naxis=2) if 'CTYPE1' in hdul[1].header else None
            vm = cls(d['v'], d['v_err'], d['mask'].astype(bool), d['width'], d['width_err'], d['flux'],
                     d['flux_err'], d['code'], ph.get('KPCPIX'), ph.get('PSFKPC'), ph.get('LINE', ''),
                     ph.get('Z'), wcs, json.loads(ph.get('PROVEN', '[]')), json.loads(ph.get('META', '{}')))
            if 'beta_x' in d:
                vm.set_source_positions(d['beta_x'], d['beta_y'], d.get('image'))
        return vm

    # --- display ------------------------------------------------------------

    def plot(self, ax=None, quantity='v', show_masked=False, colorbar=True, **imshow_kw):
        """Image of ``quantity`` ('v', 'v_err', 'width', 'flux', 'code', ...) over the kept spaxels."""
        import matplotlib.pyplot as plt
        ax = ax or plt.gca()
        a = getattr(self, quantity).astype(float)
        if not show_masked:
            a = np.where(self.mask, a, np.nan)
        if quantity == 'v' and 'cmap' not in imshow_kw:
            lim = np.nanpercentile(np.abs(a - np.nanmedian(a)), 98) if np.isfinite(a).any() else 1.
            med = np.nanmedian(a) if np.isfinite(a).any() else 0.
            imshow_kw.update(cmap='RdBu_r', vmin=med - lim, vmax=med + lim)
        im = ax.imshow(a, **imshow_kw)
        if colorbar:
            plt.colorbar(im, ax=ax, label=quantity)
        ax.set_title(f'{self.line} {quantity}', fontsize=9)
        return im

    def __repr__(self):
        src = ', source positions' if self.has_source_positions else ''
        return f'VelocityMap({self.line}, shape={self.shape}, kept={int(self.mask.sum())}{src})'


def _jd(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, tuple):
        return list(o)
    return str(o)
