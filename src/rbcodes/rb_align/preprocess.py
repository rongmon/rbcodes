"""
preprocess.py — thin wrapper around GUIs.ifuviewer.processing.cube_collapse.

Collapses IFU cubes to 2D images for alignment.
No-op for 2D images.
"""
from __future__ import annotations

from typing import Callable, List, Optional

import numpy as np

from rbcodes.GUIs.ifuviewer.processing.cube_collapse import (
    build_whitelight,
    build_narrowband,
)

from .core import Frame


def _preprocess_frame(frame: Frame,
                      method: str = 'whitelight',
                      wl_range=None,
                      func: Optional[Callable] = None) -> None:
    """
    Collapse one Frame's 3D cube to a 2D image and store as frame.image2d.

    Parameters
    ----------
    frame : Frame
    method : str
        'whitelight' | 'narrowband'
    wl_range : [wmin, wmax] or None
        Wavelength range in Angstroms; required for 'narrowband'.
    func : callable or None
        Custom collapse: func(flux_3d) -> image_2d. Bypasses built-ins.
    """
    # 2D images — no-op; image2d is already the data
    if frame.input_type == 'image' or frame.data.ndim == 2:
        if frame.image2d is None:
            frame.image2d = frame.data
        return

    flux = frame.data   # shape (n_wave, ny, nx)
    wave = frame.wave

    if func is not None:
        # User-supplied collapse function
        result = func(flux)
    elif method == 'narrowband':
        if wl_range is None:
            raise ValueError(
                "preprocess(method='narrowband') requires wl_range=[wmin, wmax]."
            )
        wmin, wmax = float(wl_range[0]), float(wl_range[1])
        result = build_narrowband(flux, wave, wmin, wmax)
    else:
        # Default: whitelight
        result = build_whitelight(flux, wave)

    frame.image2d = np.asarray(result, dtype=np.float64)
