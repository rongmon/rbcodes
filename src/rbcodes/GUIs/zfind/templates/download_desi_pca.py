"""
download_desi_pca.py
--------------------
Download DESI redrock PCA eigenvector templates from GitHub and save them
in the pca/ subdirectory alongside this script.

The DESI redrock-templates package is BSD-licensed:
  https://github.com/desihub/redrock-templates

Each FITS file contains a 2-D array of shape (n_components, n_pixels)
on a log10-linear wavelength grid defined by CRVAL1/CDELT1 in the header.

Usage
-----
    python download_desi_pca.py

Run once; re-run to refresh.  The engine loads from the saved FITS files.
"""

import os
import sys
import urllib.request

DESI_BASE = (
    "https://raw.githubusercontent.com/desihub/redrock-templates/main/"
)

WANTED = {
    'galaxy':  'rrtemplate-GALAXY-None-v2.6.fits',
    'qso_loz': 'rrtemplate-QSO-LOZ-v1.1.fits',
    'qso_hiz': 'rrtemplate-QSO-HIZ-v1.1.fits',
}

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pca')


def _fetch(url: str, dest: str) -> None:
    print(f"  Downloading {os.path.basename(dest)} …")
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = r.read()
        with open(dest, 'wb') as f:
            f.write(data)
        print(f"    Saved → {dest}  ({len(data) / 1024:.0f} kB)")
    except Exception as e:
        print(f"    ERROR fetching {url}: {e}", file=sys.stderr)
        raise


def _verify(path: str, label: str) -> None:
    """Quick sanity check: open the FITS file and print shape."""
    from astropy.io import fits
    with fits.open(path) as hdul:
        for i, hdu in enumerate(hdul):
            if hdu.data is not None and hasattr(hdu.data, 'shape'):
                hdr = hdu.header
                naxis1 = hdr.get('NAXIS1', '?')
                crval1 = hdr.get('CRVAL1', '?')
                cdelt1 = hdr.get('CDELT1', '?')
                print(f"    [{label}] HDU{i}: shape={hdu.data.shape}  "
                      f"CRVAL1={crval1}  CDELT1={cdelt1}  NAXIS1={naxis1}")
                break


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    n_ok = 0
    for label, fname in WANTED.items():
        url  = DESI_BASE + fname
        dest = os.path.join(OUTDIR, fname)
        print(f"\n{label}: {fname}")
        try:
            _fetch(url, dest)
            _verify(dest, label)
            n_ok += 1
        except Exception:
            pass

    if n_ok == 0:
        print("\nERROR: no templates downloaded.", file=sys.stderr)
        sys.exit(1)

    print(f"\nDone — {n_ok}/{len(WANTED)} PCA template sets saved to {OUTDIR}/")


if __name__ == '__main__':
    main()
