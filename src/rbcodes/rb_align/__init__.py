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
  c = wcs_align.from_file(reference='hst_f814w.fits',
                           targets=['kcwi_obs.fits'],
                           input_type='ifu')
  c.preprocess(method='whitelight', collapse='mean')
  c.find_sources(strategy='interactive', stretch='zscale', box=0.1,
                 save_catalog='sources.fits')   # save for batch reuse
  c.align()
  c.qa(plot=True)
  c.write_output()                    # → kcwi_obs_wcsfix.fits
  c.apply_to('kcwi_obs_var.fits')    # → kcwi_obs_var_wcsfix.fits

Interactive window controls
  Left-click reference   — pick a source (recentroids automatically)
  Left-click target      — confirm position for current pending source
  Double-click any marker — toggle edit mode (cyan circle appears on target;
                            left-click near the circle to re-place that source)
  Right-click            — IDLE: delete nearest pair  |  PENDING: cancel
  u                      — delete nearest pair (hover cursor near it first)
  Space                  — auto-accept predicted target position
  Enter                  — finish and close

Example 2 — Batch mode: load saved catalog, inspect, edit bad pairs
-------------------------------------------------------------------
  # Step 1: run interactively once and save the catalog
  c = wcs_align.from_file(reference='hst_f814w.fits',
                           targets=['kcwi_01.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='interactive', save_catalog='sources.fits')
  c.align()
  c.write_output()

  # Step 2: load the catalog for a new target
  #   → All pairs shown as numbered red markers.
  #   → Double-click any marker to toggle edit mode (cyan circle).
  #   → Left-click near the circle to re-place that source.
  #   → Right-click or 'u' to delete a pair entirely.
  #   → Enter when satisfied.
  c2 = wcs_align.from_file(reference='hst_f814w.fits',
                            targets=['kcwi_02.fits'],
                            input_type='ifu')
  c2.preprocess()
  c2.find_sources(strategy='batch', catalog='sources.fits', box=0.1)
  # ^^ opens interactive window for inspection after loading
  c2.align()
  c2.write_output()
  c2.apply_to('kcwi_02_var.fits')

Example 3 — Cross-corr between two KCWI cubes (different FoV/pixel scale)
--------------------------------------------------------------------------
  c = wcs_align.from_file(reference='kcwi_blue.fits',
                           targets=['kcwi_red.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='cross_corr',
                 box=0.1,
                 fwhm=8.0,             # ~1" seeing at 0.147"/px
                 threshold_sigma=3.0,
                 bg_method='mad',      # robust for QSO fields
                 max_sources=50)
  print(f"Found {len(c.pairs)} pairs")
  c.align()
  c.write_output()
  c.apply_to('kcwi_red_var.fits')

Example 4 — Cross-corr then interactive refinement
---------------------------------------------------
  c = wcs_align.from_file(reference='hst_f814w.fits',
                           targets=['kcwi_obs.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='cross_corr', fwhm=8.0, bg_method='mad')
  print(f"Initial pairs: {c.pairs}")
  c.remove_pair(1)            # drop a bad pair by index (0-based)
  c.remove_pair(-1)           # drop last pair
  c.find_sources(strategy='interactive')   # opens with remaining pairs shown
  c.align()
  c.write_output()

Example 5 — Batch survey (10 exposures, reuse interactive catalog)
------------------------------------------------------------------
  targets = [f'kcwi_{i:02d}.fits' for i in range(10)]
  c = wcs_align.from_file(reference='hst_f814w.fits',
                           targets=targets, input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='interactive', save_catalog='field_sources.fits')
  c.align_many(on_fail='skip', min_sources=2)
  c.qa_summary()
  c.write_output()

Example 6 — Lensed arc / extended source (no point sources)
-----------------------------------------------------------
  c = wcs_align.from_file(reference='hst_arc.fits',
                           targets=['muse_arc.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='knots',
                 threshold_sigma=2.0, max_knots=10, bg_method='mad')
  c.align()
  c.write_output()

Example 7 — Narrowband collapse as reference image
---------------------------------------------------
  c = wcs_align.from_file(reference='kcwi_ref.fits',
                           targets=['kcwi_obs.fits'],
                           input_type='ifu')
  c.preprocess(method='narrowband', wl_range=[4860, 4870], collapse='median')
  c.find_sources(strategy='cross_corr', fwhm=8.0, bg_method='mad')
  c.align()
  c.write_output()

Example 8 — Apply one alignment solution to multiple files
----------------------------------------------------------
  c = wcs_align.from_file(reference='hst_f814w.fits',
                           targets=['kcwi_flux.fits'],
                           input_type='ifu')
  c.preprocess()
  c.find_sources(strategy='cross_corr', fwhm=8.0)
  c.align()
  c.write_output()                              # → kcwi_flux_wcsfix.fits
  c.apply_to('kcwi_var.fits')                  # → kcwi_var_wcsfix.fits
  c.apply_to('kcwi_mask.fits',
             output='kcwi_mask_corr.fits')     # custom output name
  c.apply_to(['kcwi_var.fits', 'kcwi_mask.fits', 'kcwi_sky.fits'])

  # Access corrected WCS directly
  wcs = c.corrected_wcs          # astropy WCS for target 0
  wcs = c.get_corrected_wcs(0)   # same, explicit index

Source strategies
-----------------
  interactive  — two-panel matplotlib window; click pairs manually (always works)
                 existing pairs shown on open (after batch or remove_pair)
  cross_corr   — Hough catalog matching in sky space (works across pixel scales)
                 kwargs: fwhm, threshold_sigma, bg_method, max_sources
  dao          — blind DAOStarFinder detection (requires photutils)
                 kwargs: fwhm, threshold_sigma, bg_method, max_sources
  knots        — segmentation centroids for arcs / extended sources
                 kwargs: threshold_sigma, max_knots, bg_method
  gaia         — Gaia DR3 catalog matching (requires astroquery, wide fields)
                 kwargs: radius_deg, max_stars
  batch        — load saved FITS catalog, reproject onto current frames,
                 recentroid, then open interactive window for inspection/editing
                 kwargs: catalog (required), box
  auto         — cascade: gaia → dao → knots → cross_corr → interactive

find_sources() keyword options
-------------------------------
  box             — centroid half-width in arcsec (default 0.1"); converted to
                    pixels per-frame via WCS pixel scale (min 3 px enforced)
                    Typical: 0.05–0.5" for HST (~0.05"/px),
                             0.3–1.0" for KCWI/MUSE (~0.15–0.2"/px)
  fwhm            — source FWHM in px; 0 = auto (cross_corr, dao)
  threshold_sigma — detection threshold in units of background σ (default 3.0)
  bg_method       — 'sigma_clip' (default) | 'mad' (more robust for IFU/QSO fields)
  max_sources     — max sources to use from each image (cross_corr, dao)
  max_knots       — max knots to detect (knots strategy)
  radius_deg      — Gaia query radius in degrees (gaia strategy)
  max_stars       — max Gaia stars to use (gaia strategy)
  save_catalog    — path to save matched pairs as FITS (any strategy)
  catalog         — path to load pairs from (batch strategy, required)
  stretch         — display stretch: 'zscale' | '99%' | '95%' | 'minmax' | (vmin,vmax)

preprocess() options
---------------------
  method   — 'whitelight' (default) | 'narrowband'
  wl_range — [wmin, wmax] in Angstroms (narrowband only)
  collapse — 'mean' (default) | 'median' | 'sum'
  func     — custom collapse function: func(flux_3d) → image_2d

Pair management
---------------
  c.pairs              — list of (ref_x, ref_y, ra, dec, tgt_x, tgt_y) tuples
  c.remove_pair(i)     — remove pair by index (0-based; -1 = last)
  c.remove_last_pair() — remove last pair (undo)
  c.clear_pairs()      — remove all pairs

WCS output
----------
  c.corrected_wcs          — corrected 2D WCS for target 0 (after align())
  c.get_corrected_wcs(i)   — corrected WCS for target i
  c.apply_to(path)         — apply correction to another FITS file
  c.apply_to([p1, p2])     — apply to multiple files
  c.apply_to(p, output=q)  — custom output path
  c.write_output()         — write all corrected targets (adds _wcsfix suffix)

Full example scripts
--------------------
  See rbcodes/rb_align/examples/align_examples.py
"""
    print(msg)
