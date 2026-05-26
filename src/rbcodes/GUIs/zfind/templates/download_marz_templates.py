"""
download_marz_templates.py
--------------------------
Fetch the MARZ spectral templates from GitHub and save them as FITS files
in the marz/ subdirectory alongside this script.

MARZ is MIT-licensed: https://github.com/Samreay/Marz
Templates are stored as JavaScript arrays in templates.js; this script
parses that file and converts the chosen galaxy/QSO templates to FITS.

Usage
-----
    python download_marz_templates.py

Run once; re-run to refresh.  The engine will load from the saved FITS files.
"""

import os
import re
import sys
import urllib.request

import numpy as np
from astropy.io import fits

MARZ_JS_URL = (
    "https://raw.githubusercontent.com/Samreay/Marz/master/js/templates.js"
)

# Templates we want:  JS name → output filename
WANTED = {
    "Early Type Absorption Galaxy":  "EarlyType.fits",
    "Intermediate Type Galaxy":      "Intermediate.fits",
    "Late Type Emission Galaxy":     "LateTypeEmission.fits",
    "Composite Galaxy":              "Composite.fits",
    "Quasar":                        "QSO.fits",
    "High Redshift Star Forming Galaxy": "HighZSFG.fits",
}

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "marz")


def _fetch_js() -> str:
    print(f"Downloading {MARZ_JS_URL} …")
    with urllib.request.urlopen(MARZ_JS_URL, timeout=30) as r:
        return r.read().decode("utf-8")


def _parse_templates(js: str) -> dict:
    """Return dict: name → {'wave': ndarray, 'flux': ndarray}."""
    pattern = re.compile(
        r"name:\s*'(?P<name>[^']+)'.*?"
        r"start_lambda:\s*(?P<sl>[\d.e+\-]+).*?"
        r"end_lambda:\s*(?P<el>[\d.e+\-]+).*?"
        r"log_linear:\s*(?P<ll>true|false).*?"
        r"spec:\s*\[",
        re.DOTALL,
    )
    out = {}
    for m in pattern.finditer(js):
        name = m.group("name")
        if name not in WANTED:
            continue
        sl = float(m.group("sl"))
        el = float(m.group("el"))
        log_linear = m.group("ll") == "true"
        # read the spec array
        start = m.end()
        end = js.index("]", start)
        flux = np.array([float(x) for x in js[start:end].split(",")])
        n = len(flux)
        if log_linear:
            # start/end are log10(Angstrom)
            wave = 10.0 ** np.linspace(sl, el, n)
        else:
            # start/end are Angstrom (or may already be linear)
            # Some entries have start in log10 space even when log_linear=false;
            # detect by value: if sl < 10 it's likely log10
            if sl < 10.0:
                wave = 10.0 ** np.linspace(sl, el, n)
            else:
                wave = np.linspace(sl, el, n)
        out[name] = {"wave": wave, "flux": flux}
        print(f"  Parsed '{name}': {n} pixels, "
              f"{wave.min():.0f}–{wave.max():.0f} Å")
    return out


def _save_fits(name: str, wave: np.ndarray, flux: np.ndarray,
               filename: str) -> None:
    hdr = fits.Header()
    hdr["SIMPLE"]  = True
    hdr["TNAME"]   = name
    hdr["CRVAL1"]  = float(wave[0])
    hdr["CDELT1"]  = float(wave[1] - wave[0])   # approx; wave stored explicitly
    hdr["CTYPE1"]  = "WAVE"
    hdr["CUNIT1"]  = "Angstrom"
    hdr["SOURCE"]  = "MARZ (MIT licence) github.com/Samreay/Marz"
    primary = fits.PrimaryHDU(header=hdr)
    wave_hdu = fits.ImageHDU(data=wave.astype(np.float64), name="WAVE")
    flux_hdu = fits.ImageHDU(data=flux.astype(np.float64), name="FLUX")
    fits.HDUList([primary, wave_hdu, flux_hdu]).writeto(filename, overwrite=True)
    print(f"  Saved → {os.path.basename(filename)}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    js = _fetch_js()
    templates = _parse_templates(js)

    if not templates:
        print("ERROR: no templates parsed — check the JS URL or regex.")
        sys.exit(1)

    for js_name, fits_name in WANTED.items():
        if js_name not in templates:
            print(f"  WARNING: '{js_name}' not found in JS — skipping.")
            continue
        t = templates[js_name]
        outpath = os.path.join(OUTDIR, fits_name)
        _save_fits(js_name, t["wave"], t["flux"], outpath)

    print(f"\nDone — {len(templates)} templates saved to {OUTDIR}/")


if __name__ == "__main__":
    main()
