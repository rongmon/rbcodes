"""
io.py — thin wrappers around GUIs.ifuviewer.io for rb_align.

Responsibilities
----------------
- _load_frame_from_file  : build a Frame from a FITS path
- _frame_from_tuple      : build a Frame from (array, wcs, header)
- _write_output          : write corrected FITS files
- _update_header_inplace : patch WCS keywords in original file
- _output_path           : filename convention helper
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.io import fits

from rbcodes.GUIs.ifuviewer.io.auto_cube import load_fits
from rbcodes.GUIs.ifuviewer.io.image2d import load_image

from .core import Frame


# ---------------------------------------------------------------------------
# Frame construction
# ---------------------------------------------------------------------------

def _auto_detect_type(path: str) -> str:
    """
    Detect whether a FITS file is a 2D image or IFU cube.

    Uses NAXIS to decide: if any extension has NAXIS == 3 → 'ifu', else 'image'.
    """
    p = Path(path)
    with fits.open(p) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim == 3:
                return 'ifu'
    return 'image'


def _load_frame_from_file(path: str,
                           input_type: Optional[str] = None,
                           auto_detect: bool = False) -> Frame:
    """
    Load a FITS file and return a Frame.

    Parameters
    ----------
    path : str
        Path to the FITS file.
    input_type : str or None
        'ifu' or 'image'.  If None and auto_detect=True, NAXIS is checked.
    auto_detect : bool
        If True, ignore input_type and detect from file.
    """
    if auto_detect or input_type is None:
        detected = _auto_detect_type(path)
        input_type = detected

    if input_type == 'image':
        obj = load_image(path)
    else:
        obj = load_fits(path)
        # load_fits may return a FITSImage if no 3D data found
        actual_type = 'image' if obj.flux.ndim == 2 else 'ifu'
        input_type = actual_type

    frame = Frame(
        data=obj.flux,
        wcs=obj.wcs,           # already 2D spatial WCS from load_fits/_make_2d_wcs
        header=obj.header,
        input_type=input_type,
        wave=obj.wave,
        image2d=None,
        path=str(path),
    )
    # For 2D images, image2d is just the data itself
    if input_type == 'image':
        frame.image2d = obj.flux
    return frame


def _frame_from_tuple(tup, input_type: str) -> Frame:
    """
    Build a Frame from a (array, wcs, header) tuple.

    Parameters
    ----------
    tup : (array, wcs, header)
    input_type : str
        'ifu' or 'image'
    """
    arr, wcs, header = tup
    arr = np.asarray(arr, dtype=np.float64)
    actual_type = 'image' if arr.ndim == 2 else input_type
    frame = Frame(
        data=arr,
        wcs=wcs,
        header=header,
        input_type=actual_type,
        wave=None,
        image2d=None,
        path=None,
    )
    if actual_type == 'image':
        frame.image2d = arr
    return frame


# ---------------------------------------------------------------------------
# Output path helper
# ---------------------------------------------------------------------------

def _output_path(input_path: str, suffix: str = '_wcsfix') -> Path:
    """
    Build output file path by inserting suffix before the extension.

    Example
    -------
    _output_path('/data/kcwi_01.fits', '_wcsfix')
    → PosixPath('/data/kcwi_01_wcsfix.fits')
    """
    p = Path(input_path)
    return p.parent / f"{p.stem}{suffix}{p.suffix}"


# ---------------------------------------------------------------------------
# Writing corrected FITS files
# ---------------------------------------------------------------------------

def _update_3d_header(header_3d, new_2d_wcs) -> None:
    """
    Patch only spatial WCS keywords back into a 3D cube header in-place.

    Wavelength axis keywords (CRVAL3, CDELT3, CTYPE3 …) are never touched.
    """
    new_hdr = new_2d_wcs.to_header()
    spatial_keys = [
        'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
        'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
        'CDELT1', 'CDELT2',
        'CTYPE1', 'CTYPE2',
        'CUNIT1', 'CUNIT2',
    ]
    for key in spatial_keys:
        if key in new_hdr:
            header_3d[key] = new_hdr[key]


def _write_output(aligner, suffix: str = '_wcsfix') -> None:
    """
    Write WCS-corrected FITS files for all targets.

    For each target, the new WCS is patched into a copy of the original header,
    and the file is written with the given suffix.  Originals are untouched.
    """
    for idx, (frame, path) in enumerate(
            zip(aligner.targets, aligner._target_paths)):
        if path is None:
            print(f"[rb_align] target {idx}: no file path — skipping write.")
            continue

        new_wcs = (aligner._new_wcs[idx]
                   if idx < len(aligner._new_wcs) else None)
        if new_wcs is None:
            print(f"[rb_align] target {idx}: no corrected WCS — skipping write.")
            continue

        out_path = _output_path(path, suffix)

        # Build corrected header
        hdr = copy.deepcopy(frame.header)
        if frame.input_type == 'ifu':
            _update_3d_header(hdr, new_wcs)
        else:
            # 2D image — write full 2D WCS header
            new_hdr = new_wcs.to_header()
            for key, val in new_hdr.items():
                hdr[key] = val

        # Re-open original to preserve all extensions
        if Path(path).exists():
            with fits.open(path) as hdul:
                hdul_out = copy.deepcopy(hdul)
                # Find the extension that holds the primary data and patch header
                for ext in hdul_out:
                    if ext.data is not None and ext.data.shape == frame.data.shape:
                        for key, val in hdr.items():
                            try:
                                ext.header[key] = val
                            except Exception:
                                pass
                        break
                else:
                    # Fallback: patch primary header
                    for key, val in hdr.items():
                        try:
                            hdul_out[0].header[key] = val
                        except Exception:
                            pass
                hdul_out.writeto(str(out_path), overwrite=True)
        else:
            # No original file — write from data
            primary_hdu = fits.PrimaryHDU(data=frame.data, header=hdr)
            hdul_out = fits.HDUList([primary_hdu])
            hdul_out.writeto(str(out_path), overwrite=True)

        print(f"[rb_align] wrote {out_path}")


def _update_header_inplace(aligner, target_idx: int) -> None:
    """
    Patch WCS keywords in the original FITS file for one target in-place.

    Uses fits.update() for minimal modification.
    """
    frame = aligner.targets[target_idx]
    path = aligner._target_paths[target_idx]
    if path is None:
        raise ValueError(f"target {target_idx} has no file path.")

    new_wcs = aligner._new_wcs[target_idx] if target_idx < len(aligner._new_wcs) else None
    if new_wcs is None:
        raise ValueError(f"target {target_idx} has no corrected WCS — run align() first.")

    hdr_patch = copy.deepcopy(frame.header)
    if frame.input_type == 'ifu':
        _update_3d_header(hdr_patch, new_wcs)
    else:
        new_hdr = new_wcs.to_header()
        for key, val in new_hdr.items():
            hdr_patch[key] = val

    with fits.open(path, mode='update') as hdul:
        for ext in hdul:
            if ext.data is not None and ext.data.shape == frame.data.shape:
                for key, val in hdr_patch.items():
                    try:
                        ext.header[key] = val
                    except Exception:
                        pass
                break
        else:
            for key, val in hdr_patch.items():
                try:
                    hdul[0].header[key] = val
                except Exception:
                    pass
        hdul.flush()
    print(f"[rb_align] updated header in {path}")


# ---------------------------------------------------------------------------
# Apply corrected WCS to arbitrary external files
# ---------------------------------------------------------------------------

def _apply_wcs_to_files(new_wcs, paths: list, outputs: list,
                         suffix: str = '_wcsfix') -> None:
    """
    Patch the spatial WCS from *new_wcs* into each file in *paths* and write
    to the corresponding entry in *outputs* (None → auto-name with *suffix*).

    The data arrays are never modified — only header keywords are patched.
    All extensions are preserved.  Auto-detects 2D vs 3D from NAXIS.
    """
    for path, out in zip(paths, outputs):
        p = Path(path)
        if not p.exists():
            print(f"[rb_align] apply_to: '{path}' not found — skipping.")
            continue

        out_path = Path(out) if out is not None else _output_path(path, suffix)

        with fits.open(p) as hdul:
            hdul_out = copy.deepcopy(hdul)

        # Patch every extension that has data (covers multi-extension FITS)
        patched = False
        for ext in hdul_out:
            if ext.data is None:
                continue
            ndim = ext.data.ndim
            if ndim == 3:
                _update_3d_header(ext.header, new_wcs)
                patched = True
            elif ndim == 2:
                for key, val in new_wcs.to_header().items():
                    try:
                        ext.header[key] = val
                    except Exception:
                        pass
                patched = True

        if not patched:
            print(f"[rb_align] apply_to: no data extensions found in '{path}' — skipping.")
            continue

        hdul_out.writeto(str(out_path), overwrite=True)
        print(f"[rb_align] apply_to: wrote {out_path}")
