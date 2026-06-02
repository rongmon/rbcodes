# GUI Tools

[Back to Main Page](../main_readme.md)

## Quick Reference

| Tool | Description | Docs |
|------|-------------|------|
| **launch_specgui** | Interactive `rb_spec` GUI — single spectrum and batch modes | [rb_spec docs](rb_spec/rb_spec.md) |
| **rb_multispec** | Multi-spectrum viewer with line ID, redshift overlay, and fitting | [multispec docs](multispec/multispec.md) |
| **rb_ifuview** | IFU datacube viewer — KCWI/MUSE/generic, spectra, moment maps, ds9 | [ifuview docs](ifuview/rb_ifuview.md) |
| **rb_zgui** | Redshift measurement GUI for JWST NIRCam/grism spectroscopy | [PDF tutorial](zgui/Tutorial_for_Emission_Line_Redshift_Estimator_GUI.pdf) |
| **rb_zfind** | Semi-automated redshift finder — PCA, template, and picket-fence chi2 search; integrates with rb_zgui | [docs](zfind/rb_zfind.md) |
| **rb_llsfitter** | GUI for Lyman Limit System column density fitting | [LLSFitter docs](LLSFitter/LLSFitter.md) |
| **interactive_continuum_fit** | Polynomial and spline continuum fitter with interactive masking | [docs](interactive_continuum_fit.md) |
| **rb_align** | Astrometry alignment — IFU cubes and 2D images to a reference frame | [rb_align docs](rb_align/rb_align.md) |

---

## rb_align

[Full Documentation](rb_align/rb_align.md)

`rb_align` aligns IFU datacubes (KCWI, MUSE, NIRSpec IFU, JWST/NIRISS) and 2D
images to a reference frame using interactive or automated source matching with
full WCS fitting via `astropy`. Works standalone, in batch pipelines, and
embedded in GUIs.

```python
import rbcodes.rb_align as rb_align
rb_align.help()                         # print all workflow examples
```

```python
from rbcodes.rb_align import wcs_align

c = wcs_align.from_file(reference='hst_ref.fits', targets=['cube.fits'],
                         input_type='ifu')
c.preprocess()
c.find_sources(strategy='interactive', stretch='zscale', save_catalog='sources.fits')
c.align()
c.qa(plot=True)
c.write_output()                        # → cube_wcsfix.fits
```

**Source strategies:** `interactive` | `gaia` | `dao` | `knots` | `cross_corr` | `batch` | `auto`

**Supported modes:** image↔image, IFU↔image, IFU↔IFU, batch (10+ exposures)


