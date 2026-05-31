"""
rb_align — Astrometry alignment tool for rbcodes.

Single public entry point: wcs_align.

Usage
-----
from rbcodes.rb_align import wcs_align

c = wcs_align.from_file(reference='ref.fits', targets=['obs.fits'],
                        input_type='ifu')
c.preprocess()
c.find_sources(strategy='interactive')
c.align()
c.qa(plot=True)
c.write_output()

Call rb_align.help() for full workflow examples.
"""

from .core import wcs_align, Frame

__all__ = ['wcs_align', 'Frame', 'help']


def help():
    """Print usage examples for rb_align."""
    msg = """
rb_align — Astrometry Alignment Tool
=====================================

Quick start
-----------
  from rbcodes.rb_align import wcs_align

  c = wcs_align.from_file(reference='hst_ref.fits', targets=['cube.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='interactive')
  c.align()
  c.qa(plot=True)
  c.write_output()          # writes cube_wcsfix.fits

Example 1 — Single IFU aligned to HST (interactive)
----------------------------------------------------
  from rbcodes.rb_align import wcs_align

  c = wcs_align.from_file(reference='hst_f814w.fits',
                           targets=['kcwi_obs.fits'],
                           input_type='ifu')
  c.preprocess(method='narrowband', wl_range=[4860, 4870])
  c.find_sources(strategy='interactive', stretch='zscale', box=5,
                 save_catalog='sources.fits')   # save for batch reuse
  c.align()
  c.qa(plot=True)
  c.write_output()
  # Output: kcwi_obs_wcsfix.fits

Example 2 — Batch survey (10 exposures, reuse interactive catalog)
------------------------------------------------------------------
  from rbcodes.rb_align import wcs_align

  targets = [f'kcwi_{i:02d}.fits' for i in range(10)]
  c = wcs_align.from_file(reference='hst_f814w.fits',
                           targets=targets, input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='interactive', save_catalog='field_sources.fits')
  c.align_many(on_fail='skip', min_sources=2)
  c.qa_summary()
  c.write_output()
  # Output: kcwi_00_wcsfix.fits ... kcwi_09_wcsfix.fits

Example 3 — KCWI blue/red cross-instrument (cross-correlation)
--------------------------------------------------------------
  c = wcs_align.from_file(reference='kcwi_blue.fits',
                           targets=['kcwi_red.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='cross_corr')   # same PSF, same morphology
  c.align()
  c.write_output()
  # Output: kcwi_red_wcsfix.fits

Example 4 — Lensed arc / extended source (no point sources)
-----------------------------------------------------------
  c = wcs_align.from_file(reference='hst_arc.fits',
                           targets=['muse_arc.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='knots')   # segmentation → bright knot centroids
  c.align()
  c.qa(plot=True)
  c.write_output()

Example 5 — Load saved catalog to skip interactive step
-------------------------------------------------------
  c = wcs_align.from_file(reference='hst_f814w.fits',
                           targets=['kcwi_obs.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='batch', catalog='sources.fits')
  c.align()
  c.write_output()

Source strategies
-----------------
  interactive  — click pairs in two-panel matplotlib window (always works)
  gaia         — Gaia DR3 catalog matching (requires astroquery, wide fields)
  dao          — blind DAOStarFinder detection (requires photutils)
  knots        — segmentation centroids for arcs/extended sources
  cross_corr   — cross-correlation, same-instrument IFU↔IFU only
  batch        — reuse saved catalog (no user input)
  auto         — cascade: gaia → dao → knots → cross_corr → interactive

display stretch options for find_sources()
------------------------------------------
  'zscale' (default) | '99.5%' | '99%' | '98%' | '97%' | '95%' | 'minmax'
  or a manual tuple: stretch=(vmin, vmax)

Full example scripts
--------------------
  See rbcodes/rb_align/examples/
"""
    print(msg)
