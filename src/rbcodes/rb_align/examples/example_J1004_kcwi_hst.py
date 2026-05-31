"""
Example: Align KCWI blue cube to KCWI red cube for J1004+4112.

The red cube acts as the astrometric reference; the blue cube is the target.
Adjust the file paths below to match your local data layout.

Usage:
    python example_J1004_kcwi_hst.py
"""

from rbcodes.rb_align import wcs_align

reference_image = (
    '/Users/bordoloi/Dropbox/Research/KCWI/Lensed_QSO/J1004+4112/'
    'new-data/2024/new_reduction/J1004p4112_red/'
    'J1004p4112_red_2025_03_03_N169_Large_RL_7200s_kskywizard_zap_icubes_flux.fits'
)

ifu_cube = (
    '/Users/bordoloi/Dropbox/Research/KCWI/Lensed_QSO/J1004+4112/'
    'new-data/2024/new_reduction/J1004p4112_blue/'
    'J1004p4112_blue_2025_03_03_N169_Large_BL_7920s_icubes_flux.fits'
)

# ---------------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------------
print('Loading files...')
c = wcs_align.from_file(
    reference=reference_image,
    targets=[ifu_cube],
    input_type='ifu',
)
print(f'  Reference : {c.reference.path}')
print(f'  Target    : {c.targets[0].path}')

# ---------------------------------------------------------------------------
# 2. Collapse IFU cube to 2D white-light image for alignment.
#    Use narrowband instead if a bright line dominates, e.g.:
#      c.preprocess(method='narrowband', wl_range=[7000, 7100])
# ---------------------------------------------------------------------------
print('Preprocessing...')
c.preprocess()
print('  Done.')

# ---------------------------------------------------------------------------
# 3. Interactive source selection
#
#    A two-panel window opens (reference left, target right).
#    - Left-click on reference to pick a source.
#    - Left-click near the cyan prediction circle on the target to confirm.
#    - Space  : auto-accept prediction (good when WCS is already close).
#    - Right-click : cancel pending reference click.
#    - u      : undo last pair.
#    - Enter  : done.
#
#    Aim for 4–6 sources spread across the field for a good affine fit.
#    The catalog is saved so you can skip this step next time (see 3b).
# ---------------------------------------------------------------------------
print('Opening interactive alignment window...')
c.find_sources(
    strategy='interactive',
    stretch='zscale',               # 'zscale' | '99%' | 'minmax' | (vmin, vmax)
    box=5,                          # centroid refinement half-width in pixels
    save_catalog='J1004_sources.fits',
)
print(f'  {len(c.pairs)} source pairs selected.')

# -- 3b. Skip interactive next time by loading the saved catalog:
# c.find_sources(strategy='batch', catalog='J1004_sources.fits')

# ---------------------------------------------------------------------------
# 4. Align
# ---------------------------------------------------------------------------
print('Fitting WCS...')
c.align()
print(f'  Fit type  : {c.fit_type[0]}')
print(f'  RMS       : {c.rms_residuals[0]:.4f} arcsec')

# ---------------------------------------------------------------------------
# 5. QA
# ---------------------------------------------------------------------------
c.qa(plot=True)

# ---------------------------------------------------------------------------
# 6. Write corrected cube (same directory as input, suffix _wcsfix)
# ---------------------------------------------------------------------------
print('Writing output...')
c.write_output(suffix='_wcsfix')
print('Done.')
