"""
Generic rb_align usage examples.

Replace the file paths below with your own data.
Run rb_align.help() in Python for a compact summary of all workflows.

Usage:
    python example_generic.py
"""

from rbcodes.rb_align import wcs_align

# =============================================================================
# Example 1 — Single IFU cube aligned to 2D reference image (interactive)
# =============================================================================

c = wcs_align.from_file(
    reference='hst_f814w.fits',      # 2D reference (HST, ground-based, etc.)
    targets=['kcwi_obs.fits'],       # IFU cube(s) to align
    input_type='ifu',
)

# Collapse IFU to 2D for alignment (white-light is the default)
c.preprocess()
# Or use a narrow band around an emission line:
# c.preprocess(method='narrowband', wl_range=[6550, 6570])

# Open interactive two-panel window to pick matched source pairs.
# Save catalog so you can reuse it without clicking next time.
c.find_sources(
    strategy='interactive',
    stretch='zscale',
    box=5,
    save_catalog='sources.fits',
)

c.align()
c.qa(plot=True)
c.write_output()
# → kcwi_obs_wcsfix.fits

# =============================================================================
# Example 2 — Skip interactive step by reloading a saved catalog
# =============================================================================

# c = wcs_align.from_file(reference='hst_f814w.fits',
#                          targets=['kcwi_obs.fits'], input_type='ifu')
# c.preprocess()
# c.find_sources(strategy='batch', catalog='sources.fits')
# c.align()
# c.write_output()

# =============================================================================
# Example 3 — Batch: 10 exposures, one interactive run, rest automated
# =============================================================================

# targets = [f'kcwi_{i:02d}.fits' for i in range(10)]
# c = wcs_align.from_file(reference='hst_f814w.fits',
#                          targets=targets, input_type='ifu')
# c.preprocess()
# c.find_sources(strategy='interactive', save_catalog='field_sources.fits')
# c.align_many(on_fail='skip', min_sources=2)
# c.qa_summary()
# c.write_output()

# =============================================================================
# Example 4 — Same-instrument IFU↔IFU (cross-correlation)
# =============================================================================

# c = wcs_align.from_file(reference='kcwi_blue.fits',
#                          targets=['kcwi_red.fits'], input_type='ifu')
# c.preprocess()
# c.find_sources(strategy='cross_corr')
# c.align()
# c.write_output()

# =============================================================================
# Example 5 — Lensed arc / extended sources (no point sources)
# =============================================================================

# c = wcs_align.from_file(reference='hst_arc.fits',
#                          targets=['muse_arc.fits'], input_type='ifu')
# c.preprocess()
# c.find_sources(strategy='knots')
# c.align()
# c.qa(plot=True)
# c.write_output()
