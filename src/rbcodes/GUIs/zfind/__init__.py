"""
rb_zfind — Semi-automated redshift and absorber finder.

Emission mode: galaxy/QSO redshift estimation (integrates with rb_zgui)
Absorption mode: intervening absorber finding (integrates with rb_multispec)

Usage
-----
From command line:
    rb_zfind spectrum.fits
    rb_zfind spectrum.fits -l ISM
    rb_zfind spectrum.fits -m absorption

From Python:
    from rbcodes.GUIs.zfind.main import launch_zfind
    launch_zfind('spectrum.fits')
    launch_zfind(wave, flux, error)
    launch_zfind(rb_spectrum_obj)
"""
