"""
rb_align — command-line usage examples.

Run individual examples by calling their functions, or read this file
as a reference.  All paths are placeholders — substitute your own files.
"""
from rbcodes.rb_align import wcs_align


# ---------------------------------------------------------------------------
# Example 1 — Single IFU aligned to HST (interactive)
# ---------------------------------------------------------------------------

def example_interactive():
    """
    Align one KCWI cube to an HST reference image interactively.

    A two-panel matplotlib window opens (reference left, target right).

    Controls:
      Left-click reference    — pick a source (auto-recentroids)
      Left-click target       — confirm position for pending source
      Double-click any marker — toggle edit mode (cyan circle on target;
                                left-click near circle to re-place)
      Right-click             — IDLE: delete nearest  |  PENDING: cancel
      u                       — delete nearest pair (hover cursor near it)
      Space                   — auto-accept predicted target position
      Enter                   — finish and close

    Aim for 4–6 sources spread across the field.
    box is in arcsec; set smaller for HST (~0.05"/px), larger for KCWI/MUSE.
    """
    c = wcs_align.from_file(
        reference='hst_f814w.fits',
        targets=['kcwi_obs.fits'],
        input_type='ifu')
    c.preprocess(method='whitelight')
    c.find_sources(strategy='interactive', stretch='zscale', box=0.1,
                   save_catalog='sources.fits')   # save for batch reuse
    c.align()
    c.qa(plot=True)
    c.write_output()                    # → kcwi_obs_wcsfix.fits

    # Apply the same correction to the variance cube — no re-alignment
    c.apply_to('kcwi_obs_var.fits')    # → kcwi_obs_var_wcsfix.fits


# ---------------------------------------------------------------------------
# Example 2 — Batch mode: load catalog, inspect, edit bad pairs
# ---------------------------------------------------------------------------

def example_batch():
    """
    Reload a previously saved catalog onto a new target.

    After loading, the interactive window opens automatically with all
    catalog pairs shown as numbered red markers.  Pairs whose positions
    look wrong can be corrected without re-doing the whole session:

      Double-click a marker   — enter edit mode (cyan circle on target)
      Left-click near circle  — re-place that source (recentroids)
      Right-click / u         — delete the pair entirely
      Enter                   — accept and close

    box controls the centroid refinement radius when reprojecting the
    catalog positions onto the current reference and target frames.
    Use a larger value for KCWI/MUSE (0.3–1.0") vs HST (0.05–0.3").
    """
    # Session 1: interactive alignment, save catalog
    c = wcs_align.from_file(
        reference='hst_f814w.fits',
        targets=['kcwi_01.fits'],
        input_type='ifu')
    c.preprocess()
    c.find_sources(strategy='interactive', save_catalog='sources.fits')
    c.align()
    c.write_output()
    c.apply_to('kcwi_01_var.fits')

    # Session 2: load catalog for a new target, opens interactive for inspection
    c2 = wcs_align.from_file(
        reference='hst_f814w.fits',
        targets=['kcwi_02.fits'],
        input_type='ifu')
    c2.preprocess()
    c2.find_sources(strategy='batch', catalog='sources.fits', box=0.1)
    # ^^ interactive window opens automatically — inspect and edit if needed
    c2.align()
    c2.write_output()
    c2.apply_to('kcwi_02_var.fits')


# ---------------------------------------------------------------------------
# Example 3 — Cross-correlation (Hough catalog, multi-source)
# ---------------------------------------------------------------------------

def example_cross_corr():
    """
    Align two KCWI cubes with different plate scales / FoVs.
    Detects sources independently in each image, votes for the best shift
    via a 2-D Hough histogram, and returns one pair per matched source.

    Options tunable via find_sources kwargs:
      fwhm            — source FWHM in px (0 = auto)
      threshold_sigma — detection threshold in units of background σ
      bg_method       — 'sigma_clip' (default) or 'mad' (more robust for
                        IFU data with correlated noise or high source fraction)
      max_sources     — max sources to use from each image
    """
    c = wcs_align.from_file(
        reference='kcwi_blue.fits',
        targets=['kcwi_red.fits'],
        input_type='ifu')
    c.preprocess()
    c.find_sources(strategy='cross_corr',
                   box=0.1,
                   fwhm=8.0,              # ~1" seeing at 0.147"/px
                   threshold_sigma=3.0,
                   bg_method='mad',       # robust for QSO fields
                   max_sources=50)
    print(f"Found {len(c.pairs)} pairs")
    c.align()
    c.qa(plot=True)
    c.write_output()
    c.apply_to('kcwi_red_var.fits')


# ---------------------------------------------------------------------------
# Example 4 — Cross-corr then interactive refinement
# ---------------------------------------------------------------------------

def example_cross_corr_then_interactive():
    """
    Use cross_corr to get initial pairs, then open the interactive window
    to verify, add, or remove pairs before aligning.

    The interactive window opens with all cross_corr pairs already shown
    as numbered red markers.  Double-click any marker to toggle edit mode
    (cyan circle on target), then left-click to re-place it.
    """
    c = wcs_align.from_file(
        reference='hst_f814w.fits',
        targets=['kcwi_obs.fits'],
        input_type='ifu')
    c.preprocess()

    # Step 1: automated initial guess
    c.find_sources(strategy='cross_corr', fwhm=8.0, bg_method='mad')
    print(f"cross_corr found {len(c.pairs)} pairs: {c.pairs}")

    # Step 2: programmatically drop a known-bad pair if needed
    # c.remove_pair(1)       # remove pair #2 (0-based)

    # Step 3: open interactive — existing pairs shown, edit / add / delete freely
    c.find_sources(strategy='interactive')

    c.align()
    c.write_output()
    c.apply_to(['kcwi_obs_var.fits', 'kcwi_obs_mask.fits'])


# ---------------------------------------------------------------------------
# Example 5 — DAO point-source detection
# ---------------------------------------------------------------------------

def example_dao():
    """
    Blind point-source detection with DAOStarFinder (photutils) or fallback.
    Detects sources in the reference, projects via WCS to the target.
    Best for fields with multiple well-resolved point sources.
    """
    c = wcs_align.from_file(
        reference='hst_f814w.fits',
        targets=['kcwi_obs.fits'],
        input_type='ifu')
    c.preprocess()
    c.find_sources(strategy='dao',
                   box=0.1,
                   fwhm=6.0,
                   threshold_sigma=3.0,
                   bg_method='sigma_clip',
                   max_sources=30)
    c.align()
    c.write_output()


# ---------------------------------------------------------------------------
# Example 6 — Extended source / lensed arc (knots)
# ---------------------------------------------------------------------------

def example_knots():
    """
    Segmentation-based knot detection for arcs or extended sources
    with no point sources.
    """
    c = wcs_align.from_file(
        reference='hst_arc.fits',
        targets=['muse_arc.fits'],
        input_type='ifu')
    c.preprocess()
    c.find_sources(strategy='knots',
                   box=0.1,
                   threshold_sigma=2.0,
                   max_knots=10,
                   bg_method='mad')
    c.align()
    c.write_output()


# ---------------------------------------------------------------------------
# Example 7 — Narrowband / median collapse for reference
# ---------------------------------------------------------------------------

def example_narrowband_ref():
    """
    Collapse the reference cube in a specific wavelength range before aligning.
    Useful when a bright emission line is the best astrometric signal.
    """
    c = wcs_align.from_file(
        reference='kcwi_ref.fits',
        targets=['kcwi_obs.fits'],
        input_type='ifu')
    c.preprocess(method='narrowband', wl_range=[4860, 4870], collapse='median')
    c.find_sources(strategy='cross_corr', fwhm=8.0, bg_method='mad')
    c.align()
    c.write_output()


# ---------------------------------------------------------------------------
# Example 8 — apply_to: propagate one solution to multiple files
# ---------------------------------------------------------------------------

def example_apply_to():
    """
    After aligning the flux cube, apply the identical WCS correction to
    variance, mask, and any other co-spatial files without re-running
    source detection or alignment.
    """
    c = wcs_align.from_file(
        reference='hst_f814w.fits',
        targets=['kcwi_flux.fits'],
        input_type='ifu')
    c.preprocess()
    c.find_sources(strategy='cross_corr', fwhm=8.0, bg_method='mad')
    c.align()
    c.write_output()                     # → kcwi_flux_wcsfix.fits

    # Apply same correction — no file I/O limit, all extensions patched
    c.apply_to('kcwi_var.fits')                              # → kcwi_var_wcsfix.fits
    c.apply_to('kcwi_mask.fits', output='kcwi_mask_corr.fits')
    c.apply_to(['kcwi_var.fits', 'kcwi_mask.fits', 'kcwi_sky.fits'])

    # Access the WCS object directly for custom use
    wcs = c.corrected_wcs                # astropy WCS for target 0
    wcs2 = c.get_corrected_wcs(0)        # same, explicit index

    # Programmatic pair management before aligning
    c2 = wcs_align.from_file(
        reference='hst_f814w.fits',
        targets=['kcwi_flux.fits'],
        input_type='ifu')
    c2.preprocess()
    c2.find_sources(strategy='cross_corr', fwhm=8.0)
    print(c2.pairs)
    c2.remove_pair(1)      # drop pair index 1 (0-based)
    c2.remove_pair(-1)     # drop last pair
    c2.find_sources(strategy='interactive')   # refine remaining
    c2.align()
    c2.write_output()
