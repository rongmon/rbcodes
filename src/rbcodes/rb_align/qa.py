"""
qa.py — Quality-assurance summaries and figures for rb_align.

Public functions (called via wcs_align methods)
-----------------------------------------------
_qa(aligner, plot=False, save=None)
_qa_summary(aligner)
"""
from __future__ import annotations

from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Per-frame residuals
# ---------------------------------------------------------------------------

def _qa(aligner, plot: bool = False, save: Optional[str] = None):
    """
    Print per-frame residuals.  Optionally show or save an overlay figure.

    Parameters
    ----------
    aligner : wcs_align
    plot : bool
        If True, show the figure interactively.
    save : str or None
        Path to save the figure (PNG/PDF).  None = don't save.
    """
    n_targets = len(aligner.targets)
    if n_targets == 0:
        print("[rb_align] QA: no targets loaded.")
        return

    print("\n" + "=" * 60)
    print("  rb_align QA — per-frame residuals")
    print("=" * 60)

    for idx in range(n_targets):
        fit_type  = (aligner.fit_type[idx]
                     if idx < len(aligner.fit_type) else None)
        n_used    = (aligner.n_sources_used[idx]
                     if idx < len(aligner.n_sources_used) else 0)
        rms       = (aligner.rms_residuals[idx]
                     if idx < len(aligner.rms_residuals) else None)
        residuals = (aligner.source_residuals[idx]
                     if idx < len(aligner.source_residuals) else None)
        flagged   = idx in aligner.flagged_frames

        path_str = (aligner._target_paths[idx]
                    if idx < len(aligner._target_paths) else '<unknown>')

        status = "FLAGGED" if flagged else "OK"
        rms_str = f"{rms:.4f}\"" if rms is not None else "n/a"
        ft_str  = fit_type if fit_type is not None else "n/a"

        print(f"\n  Target {idx}: {path_str}")
        print(f"    status    : {status}")
        print(f"    fit_type  : {ft_str}")
        print(f"    n_sources : {n_used}")
        print(f"    RMS       : {rms_str}")

        if (residuals is not None and len(residuals) > 0):
            for i, r in enumerate(residuals):
                print(f"      source {i+1:>2d}: {r:.4f}\"")

    print("=" * 60 + "\n")

    if plot or save is not None:
        _make_qa_figure(aligner, show=plot, save=save)


# ---------------------------------------------------------------------------
# Compact summary table
# ---------------------------------------------------------------------------

def _qa_summary(aligner):
    """Print a compact table with one row per target."""
    n_targets = len(aligner.targets)
    if n_targets == 0:
        print("[rb_align] QA summary: no targets loaded.")
        return

    # Column widths
    col_idx  = 6
    col_path = 30
    col_fit  = 10
    col_n    = 8
    col_rms  = 12
    col_flag = 8

    rms_hdr = 'RMS(")'
    header = (
        f"{'Idx':>{col_idx}}  "
        f"{'File':<{col_path}}  "
        f"{'Fit':<{col_fit}}  "
        f"{'N_src':>{col_n}}  "
        f"{rms_hdr:{col_rms}}  "
        f"{'Flag':<{col_flag}}"
    )
    sep = "-" * len(header)
    print("\n" + sep)
    print(header)
    print(sep)

    for idx in range(n_targets):
        fit_type = (aligner.fit_type[idx]
                    if idx < len(aligner.fit_type) else None) or 'n/a'
        n_used   = (aligner.n_sources_used[idx]
                    if idx < len(aligner.n_sources_used) else 0) or 0
        rms      = (aligner.rms_residuals[idx]
                    if idx < len(aligner.rms_residuals) else None)
        flagged  = idx in aligner.flagged_frames
        path_str = (aligner._target_paths[idx]
                    if idx < len(aligner._target_paths) else '')

        if path_str is None:
            path_str = '<data>'
        # Show only the filename
        import os
        path_str = os.path.basename(str(path_str))

        rms_str  = f"{rms:.4f}" if rms is not None else "n/a"
        flag_str = "FLAGGED" if flagged else ""

        print(
            f"{idx:>{col_idx}}  "
            f"{path_str:<{col_path}}  "
            f"{fit_type:<{col_fit}}  "
            f"{n_used:>{col_n}}  "
            f"{rms_str:>{col_rms}}  "
            f"{flag_str:<{col_flag}}"
        )

    print(sep + "\n")


# ---------------------------------------------------------------------------
# QA figure
# ---------------------------------------------------------------------------

def _make_qa_figure(aligner, show: bool = True, save: Optional[str] = None):
    """
    Create a QA overlay figure.

    For each successfully aligned target, shows the target image with
    the matched source positions and residual vectors overlaid.
    """
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from .sources import _compute_norm

    n_targets = len(aligner.targets)
    n_cols = min(3, n_targets)
    n_rows = (n_targets + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(5 * n_cols, 5 * n_rows),
                              squeeze=False)

    for idx in range(n_targets):
        row, col = divmod(idx, n_cols)
        ax = axes[row][col]

        frame     = aligner.targets[idx]
        img       = frame.image2d
        residuals = (aligner.source_residuals[idx]
                     if idx < len(aligner.source_residuals) else None)
        fit_type  = (aligner.fit_type[idx]
                     if idx < len(aligner.fit_type) else 'n/a') or 'n/a'
        rms       = (aligner.rms_residuals[idx]
                     if idx < len(aligner.rms_residuals) else None)

        if img is None:
            ax.text(0.5, 0.5, f'Target {idx}\n(no image)',
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Target {idx}')
            continue

        norm = _compute_norm(img, 'zscale')
        ax.imshow(img, origin='lower', cmap='gray', norm=norm,
                  interpolation='nearest')

        # Overlay matched source positions
        pairs = aligner.pairs if idx == 0 else []
        # For batch, we don't have per-frame pairs stored separately;
        # we at least show the first-target pairs for target 0.
        for i, pair in enumerate(pairs):
            rx, ry, ra, dec, tx, ty = pair
            ax.plot(tx, ty, 'r+', ms=12, mew=2)
            ax.text(tx + 3, ty + 3, str(i + 1), color='red', fontsize=8)
            # Predicted position from current WCS (before correction)
            try:
                pred_pix = frame.wcs.all_world2pix([[ra, dec]], 0)
                px_pred, py_pred = float(pred_pix[0, 0]), float(pred_pix[0, 1])
                ax.annotate('', xy=(tx, ty), xytext=(px_pred, py_pred),
                             arrowprops=dict(arrowstyle='->', color='cyan',
                                             lw=1.5))
            except Exception:
                pass

        rms_str = f"{rms:.4f}\"" if rms is not None else "n/a"
        ax.set_title(f'Target {idx} | {fit_type} | RMS={rms_str}')
        ax.axis('off')

    # Hide unused panels
    for idx in range(n_targets, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row][col].set_visible(False)

    fig.suptitle('rb_align QA', fontsize=12)
    plt.tight_layout()

    if save is not None:
        fig.savefig(save, dpi=150, bbox_inches='tight')
        print(f"[rb_align] QA figure saved → {save}")

    if show:
        plt.show()
    else:
        plt.close(fig)
