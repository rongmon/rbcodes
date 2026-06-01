"""
core.py — wcs_align class and Frame dataclass.

The Frame dataclass holds a single image or IFU cube with its 2D spatial WCS.
wcs_align orchestrates the full pipeline: load → preprocess → find_sources
→ align → qa → write_output.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Frame dataclass
# ---------------------------------------------------------------------------

@dataclass
class Frame:
    """
    Internal container for one image or IFU cube.

    Attributes
    ----------
    data : np.ndarray
        2D image array or 3D cube (n_wave, ny, nx).
    wcs : astropy.wcs.WCS
        Always the 2D spatial WCS — never the raw 3D cube WCS.
    header : astropy.io.fits.Header
        Original FITS header (3D for IFU cubes).
    input_type : str
        'image' | 'ifu'
    wave : np.ndarray or None
        Wavelength array in Angstroms; None for 2D images.
    image2d : np.ndarray or None
        Collapsed 2D image used for alignment; filled by preprocess().
    path : str or None
        Original file path; used for write_output().
    """
    data: Any                          # np.ndarray (2D or 3D)
    wcs: Any                           # astropy.wcs.WCS, always 2D
    header: Any                        # astropy.io.fits.Header
    input_type: str                    # 'image' | 'ifu'
    wave: Optional[Any] = None         # np.ndarray or None
    image2d: Optional[Any] = None      # np.ndarray or None
    path: Optional[str] = None         # original file path

    @property
    def shape2d(self) -> Tuple[int, int]:
        """Return (ny, nx) of the 2D alignment image."""
        img = self.image2d if self.image2d is not None else (
            self.data if self.data.ndim == 2 else None
        )
        if img is None:
            raise RuntimeError("Frame has no 2D image yet — call preprocess() first.")
        return img.shape

    @property
    def ny(self) -> int:
        return self.shape2d[0]

    @property
    def nx(self) -> int:
        return self.shape2d[1]


# ---------------------------------------------------------------------------
# wcs_align class
# ---------------------------------------------------------------------------

class wcs_align:
    """
    Astrometry alignment pipeline for rbcodes.

    Typical use
    -----------
    c = wcs_align.from_file(reference='hst.fits', targets=['kcwi.fits'],
                            input_type='ifu')
    c.preprocess()
    c.find_sources(strategy='interactive')
    c.align()
    c.qa(plot=True)
    c.write_output()
    """

    # ------------------------------------------------------------------
    # Constructors (delegated to io.py)
    # ------------------------------------------------------------------

    @classmethod
    def from_file(cls,
                  reference: str,
                  targets,
                  input_type: str = 'ifu') -> 'wcs_align':
        """
        Load from FITS files.

        Parameters
        ----------
        reference : str
            Path to the reference image or cube.  Auto-detected (NAXIS check).
        targets : str or list of str
            Path(s) to target file(s).  input_type applies to all of them.
        input_type : str
            'ifu' or 'image' — applies to targets.  Reference is auto-detected.
        """
        from .io import _load_frame_from_file
        if isinstance(targets, str):
            targets = [targets]
        obj = cls.__new__(cls)
        obj._init_state()
        obj.reference = _load_frame_from_file(reference, auto_detect=True)
        obj.targets = [_load_frame_from_file(t, input_type=input_type)
                       for t in targets]
        obj._target_paths = list(targets)
        return obj

    @classmethod
    def from_data(cls,
                  reference: Tuple,
                  targets,
                  input_type: str = 'ifu') -> 'wcs_align':
        """
        Build from in-memory (array, wcs, header) tuples.

        Parameters
        ----------
        reference : (array, wcs, header)
        targets : list of (array, wcs, header)
        input_type : str
            'ifu' or 'image' — applies to all targets.
        """
        from .io import _frame_from_tuple
        if not isinstance(targets, (list, tuple)) or (
                len(targets) > 0 and not isinstance(targets[0], (list, tuple))):
            targets = [targets]
        obj = cls.__new__(cls)
        obj._init_state()
        # auto-detect reference type from data ndim
        ref_arr, ref_wcs, ref_hdr = reference
        ref_type = 'image' if np.asarray(ref_arr).ndim == 2 else 'ifu'
        obj.reference = _frame_from_tuple(reference, ref_type)
        obj.targets = [_frame_from_tuple(t, input_type) for t in targets]
        obj._target_paths = [None] * len(obj.targets)
        return obj

    # ------------------------------------------------------------------
    # Internal init helper
    # ------------------------------------------------------------------

    def _init_state(self):
        """Initialise all pipeline attributes."""
        self.reference: Frame = None
        self.targets: List[Frame] = []
        self._target_paths: List[Optional[str]] = []

        # Pair storage — list of (ref_x, ref_y, ra, dec, tgt_x, tgt_y)
        self.pairs: List[Tuple] = []

        # GUI callback
        self.on_pair_callback: Optional[Callable] = None

        # QA / align results (one entry per target)
        self.rms_residuals: List[Optional[float]] = []
        self.source_residuals: List[Optional[Any]] = []
        self.fit_type: List[Optional[str]] = []
        self.n_sources_used: List[int] = []
        self.flagged_frames: List[int] = []

        # Updated WCS objects per target (filled by align / align_many)
        self._new_wcs: List[Optional[Any]] = []

    # ------------------------------------------------------------------
    # Preprocessing (delegates to preprocess.py)
    # ------------------------------------------------------------------

    def preprocess(self,
                   method: str = 'whitelight',
                   wl_range=None,
                   func: Optional[Callable] = None):
        """
        Collapse IFU cubes to 2D images for alignment.

        Parameters
        ----------
        method : str
            'whitelight' (default) or 'narrowband'.
        wl_range : [wmin, wmax] or None
            Wavelength range in Angstroms (narrowband only).
        func : callable or None
            Custom collapse function: func(flux_3d) -> image_2d.
        """
        from .preprocess import _preprocess_frame
        # Reference
        _preprocess_frame(self.reference, method=method,
                          wl_range=wl_range, func=func)
        # Targets
        for frame in self.targets:
            _preprocess_frame(frame, method=method,
                              wl_range=wl_range, func=func)

    # ------------------------------------------------------------------
    # Source detection (delegates to sources.py)
    # ------------------------------------------------------------------

    def find_sources(self,
                     strategy: str = 'interactive',
                     stretch: str = 'zscale',
                     box: int = 5,
                     save_catalog: Optional[str] = None,
                     catalog: Optional[str] = None,
                     **kwargs):
        """
        Detect and match sources between reference and target(s).

        Parameters
        ----------
        strategy : str
            'interactive' | 'gaia' | 'dao' | 'knots' |
            'cross_corr' | 'batch' | 'auto'
        stretch : str or (vmin, vmax)
            Display stretch for interactive mode.
        box : int
            Half-width in pixels for centroid refinement.
        save_catalog : str or None
            Path to save FITS catalog after interactive session.
        catalog : str or None
            Path to existing FITS catalog (batch strategy).
        **kwargs
            Extra per-strategy options forwarded to the detection function:
            fwhm, threshold_sigma, max_sources  (cross_corr, dao)
            threshold_sigma, max_knots           (knots)
            radius_deg, max_stars                (gaia)
        """
        from .sources import _run_strategy
        _run_strategy(self, strategy=strategy, stretch=stretch,
                      box=box, save_catalog=save_catalog, catalog=catalog,
                      **kwargs)

    # ------------------------------------------------------------------
    # Alignment (delegates to align.py)
    # ------------------------------------------------------------------

    def align(self):
        """Align single target using stored pairs."""
        from .align import _align_single
        if not self.targets:
            raise RuntimeError("No targets loaded.")
        # Ensure result lists have the right length
        self._ensure_result_lists()
        _align_single(self, target_idx=0)

    def align_many(self,
                   on_fail: str = 'skip',
                   min_sources: int = 2):
        """
        Align all targets in batch.

        Parameters
        ----------
        on_fail : str
            'skip' | 'raise' | 'degrade'
        min_sources : int
            Minimum number of in-bounds sources required per frame.
        """
        from .align import _align_many
        self._ensure_result_lists()
        _align_many(self, on_fail=on_fail, min_sources=min_sources)

    def _ensure_result_lists(self):
        n = len(self.targets)
        for attr in ('rms_residuals', 'source_residuals', 'fit_type',
                     'n_sources_used', '_new_wcs'):
            lst = getattr(self, attr)
            while len(lst) < n:
                lst.append(None)

    # ------------------------------------------------------------------
    # QA (delegates to qa.py)
    # ------------------------------------------------------------------

    def qa(self,
           plot: bool = False,
           save: Optional[str] = None):
        """
        Print per-frame residuals; optionally show/save overlay figure.
        """
        from .qa import _qa
        _qa(self, plot=plot, save=save)

    def qa_summary(self):
        """Print compact table of all frames (useful after align_many)."""
        from .qa import _qa_summary
        _qa_summary(self)

    # ------------------------------------------------------------------
    # Output (delegates to io.py)
    # ------------------------------------------------------------------

    def write_output(self, suffix: str = '_wcsfix'):
        """
        Write WCS-corrected FITS files for all targets.
        Original files are never overwritten.
        """
        from .io import _write_output
        _write_output(self, suffix=suffix)

    def update_header(self, target_idx: int):
        """
        Patch the WCS keywords in the original FITS file for one target in-place.
        """
        from .io import _update_header_inplace
        _update_header_inplace(self, target_idx=target_idx)

    # ------------------------------------------------------------------
    # GUI hooks
    # ------------------------------------------------------------------

    def get_display_image(self, idx: int = 0) -> np.ndarray:
        """
        Return the 2D image for display (reference=0, or target by index).

        Returns image2d if set, else the data array (2D or collapsed 3D).
        """
        if idx == 0:
            frame = self.reference
        else:
            frame = self.targets[idx - 1]
        if frame.image2d is not None:
            return frame.image2d
        if frame.data.ndim == 2:
            return frame.data
        # Fallback: collapse cube
        return np.nanmean(frame.data, axis=0)

    def add_pair(self,
                 ref_x: float, ref_y: float,
                 tgt_x: float, tgt_y: float):
        """
        Add a matched pair by pixel coordinates.

        RA/Dec are computed from (ref_x, ref_y) via the reference WCS.
        Fires on_pair_callback if set.
        """
        wcs_ref = self.reference.wcs
        sky = wcs_ref.all_pix2world([[ref_x, ref_y]], 0)
        ra, dec = float(sky[0, 0]), float(sky[0, 1])
        self.pairs.append((ref_x, ref_y, ra, dec, tgt_x, tgt_y))
        if self.on_pair_callback is not None:
            self.on_pair_callback(ref_x, ref_y, tgt_x, tgt_y, ra, dec)

    def remove_last_pair(self):
        """Remove the last stored pair (undo)."""
        if self.pairs:
            self.pairs.pop()

    def remove_pair(self, index: int):
        """
        Remove a specific pair by index (0-based).

        Negative indices work as usual (e.g. -1 removes the last pair).

        Example
        -------
        c.find_sources(strategy='cross_corr')
        print(c.pairs)        # inspect — decide pair 1 is wrong
        c.remove_pair(1)      # remove it
        c.find_sources(strategy='interactive')   # refine the rest
        """
        if not self.pairs:
            raise IndexError("No pairs to remove.")
        try:
            del self.pairs[index]
        except IndexError:
            raise IndexError(
                f"Pair index {index} out of range (0–{len(self.pairs)-1}).")

    def clear_pairs(self):
        """Remove all stored pairs."""
        self.pairs.clear()
