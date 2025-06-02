# rbcodes/GUIs/specgui/spectrum_plotter.py
"""
Spectrum plotting module for rb_spec GUI applications.

This module provides reusable plotting functions that generate matplotlib figures
from rb_spec objects. All plots are generated fresh from the rb_spec object state
to ensure consistency and avoid update artifacts.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def plot_spectrum_overview(rb_spec, item, figure=None, clear_figure=True):
    """
    Create a comprehensive two-panel spectrum plot.
    
    Parameters
    ----------
    rb_spec : rb_spec object
        The spectrum object containing all data and analysis results
    item : object with template and analysis attributes
        Batch item containing metadata (filename, redshift, transition info, etc.)
    figure : matplotlib.figure.Figure, optional
        Existing figure to plot into. If None, creates new figure
    clear_figure : bool, optional
        Whether to clear the figure before plotting (default True)
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot
    """
    # Create or use existing figure
    if figure is None:
        fig = Figure(figsize=(10, 8))
    else:
        fig = figure
        if clear_figure:
            fig.clear()
    
    # Create subplots with shared x-axis and no space between them
    ax1, ax2 = fig.subplots(nrows=2, ncols=1, sharex=True)
    
    # Remove horizontal space between panels
    fig.subplots_adjust(hspace=0)
    
    # Top panel: Raw flux + continuum + masks (no x-axis labels)
    _plot_flux_continuum_panel(rb_spec, item, ax1, show_xlabel=False)
    
    # Bottom panel: Normalized flux + EW region + results (with x-axis labels)
    _plot_normalized_panel(rb_spec, item, ax2, show_xlabel=True)
    
    # No overall title - info is shown in status bar
    
    return fig


def _plot_flux_continuum_panel(rb_spec, item, ax, show_xlabel=False):
    """Plot the top panel: flux, error, continuum, and masked regions."""
    # Clear the axes
    ax.clear()
    
    # Plot flux and error
    ax.step(rb_spec.velo, rb_spec.flux_slice, 'k-', where='mid', 
           linewidth=0.5, label='Flux', alpha=0.8)
    
    if hasattr(rb_spec, 'error_slice'):
        ax.step(rb_spec.velo, rb_spec.error_slice, 'r-', where='mid', 
               alpha=0.6, linewidth=0.5, label='Error')
    
    # Plot continuum
    if hasattr(rb_spec, 'cont'):
        ax.plot(rb_spec.velo, rb_spec.cont, 'g-', linewidth=0.75, 
               label='Continuum', alpha=0.9)
    
    # Add continuum masks if available
    if hasattr(item, 'analysis') and hasattr(item.analysis, 'continuum_masks'):
        mask_count = 0
        for vmin, vmax in item.analysis.continuum_masks:
            label = 'Masked' if mask_count == 0 else ""
            ax.axvspan(vmin, vmax, alpha=0.2, color='red', label=label)
            mask_count += 1
    
    # Reference lines
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5, linewidth=0.5)
    ax.axvline(x=0, color='gray', linestyle=':', alpha=0.7, linewidth=0.5)
    
    # Set axis limits with buffer
    x_buffer = (max(rb_spec.velo) - min(rb_spec.velo)) * 0.05
    ax.set_xlim(min(rb_spec.velo) - x_buffer, max(rb_spec.velo) + x_buffer)
    
    # Set y limits for flux
    flux_data = rb_spec.flux_slice
    y_min = max(-0.2, min(flux_data) - 0.1 * abs(max(flux_data)))
    y_max = max(flux_data) + 0.3 * abs(max(flux_data))
    ax.set_ylim(y_min, y_max)
    
    # Add inward ticks
    ax.tick_params(direction='in', which='both')
    ax.minorticks_on()
    
    # Labels and formatting
    ax.set_ylabel('Flux', fontsize=11)
    if show_xlabel:
        ax.set_xlabel('Velocity (km/s)', fontsize=11)
    else:
        # Hide x-axis tick labels for top panel
        ax.tick_params(labelbottom=False)
    
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # Smaller legend
    #ax.legend(loc='upper right', fontsize=8)


def _plot_normalized_panel(rb_spec, item, ax, show_xlabel=True):
    """Plot the bottom panel: normalized flux, EW region, and results."""
    # Clear the axes
    ax.clear()
    
    # Plot normalized flux and error
    if hasattr(rb_spec, 'fnorm'):
        ax.step(rb_spec.velo, rb_spec.fnorm, 'k-', where='mid', 
               linewidth=0.5, label='Normalized Flux')
    
    if hasattr(rb_spec, 'enorm'):
        ax.step(rb_spec.velo, rb_spec.enorm, 'r-', where='mid', 
               alpha=0.6, linewidth=0.5, label='Error')
    
    # Reference lines
    ax.axhline(y=1.0, color='g', linestyle='--', alpha=0.7, linewidth=1)
    ax.axhline(y=0.0, color='gray', linestyle='--', alpha=0.5, linewidth=0.5)
    ax.axvline(x=0, color='gray', linestyle=':', alpha=0.7, linewidth=0.5)
    
    # EW integration region
    if hasattr(item, 'template'):
        ew_vmin = item.template.ew_vmin
        ew_vmax = item.template.ew_vmax
        
        ax.axvspan(ew_vmin, ew_vmax, alpha=0.15, color='blue', label='EW Region')
        ax.axvline(x=ew_vmin, color='b', linestyle=':', alpha=0.8, linewidth=0.8)
        ax.axvline(x=ew_vmax, color='b', linestyle=':', alpha=0.8, linewidth=0.8)
    
    # Set axis limits
    x_buffer = (max(rb_spec.velo) - min(rb_spec.velo)) * 0.05
    ax.set_xlim(min(rb_spec.velo) - x_buffer, max(rb_spec.velo) + x_buffer)
    
    # Set y limits for normalized flux
    if hasattr(rb_spec, 'fnorm'):
        norm_flux = rb_spec.fnorm
        y_min = max(-0.2, min(norm_flux) - 0.1)
        y_max = min(2.5, max(norm_flux) + 0.4)
        ax.set_ylim(y_min, y_max)
    
    # Add inward ticks
    ax.tick_params(direction='in', which='both')
    ax.minorticks_on()
    
    # Add results text box if we have measurements
    if hasattr(item, 'results') and hasattr(item.results, 'W') and item.results.W > 0:
        _add_results_text_box(item.results, ax)
    
    # Labels and formatting
    if show_xlabel:
        ax.set_xlabel('Velocity (km/s)', fontsize=11)
    ax.set_ylabel('Normalized Flux', fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # Smaller legend
    #ax.legend(loc='upper left', fontsize=8)


def _add_results_text_box(results, ax):
    """Add measurement results text box to the plot."""
    # Build results text
    ew = results.W
    ew_err = results.W_e
    logn = results.logN
    logn_err = results.logN_e
    
    info_text = f"EW = {ew:.3f} ± {ew_err:.3f} Å\n"
    

    # Handle negative N case for display
    if hasattr(results, 'N') and results.N < 0 and results.N_e > 0:
        logn = 0.0
        logn_err = np.log10(results.N_e)

    info_text += f"log N = {logn:.2f} ± {logn_err:.2f}"
    
    # Add SNR if available
    if hasattr(results, 'SNR') and results.SNR > 0:
        info_text += f"\nSNR = {results.SNR:.1f}"
    
    # Add velocity centroid and dispersion if available
    if hasattr(results, 'vel_centroid') and results.vel_centroid != 0:
        info_text += f"\nv_c = {results.vel_centroid:.1f} km/s"
    
    if hasattr(results, 'vel_disp') and results.vel_disp != 0:
        info_text += f"\nσ_v = {results.vel_disp:.1f} km/s"
    
    # Place smaller text box in bottom right
    ax.text(0.98, 0.05, info_text, 
           transform=ax.transAxes, 
           horizontalalignment='right', 
           verticalalignment='bottom',
           bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.2'),
           fontsize=8)


def plot_single_panel_normalized(rb_spec, item, figure=None, clear_figure=True):
    """
    Create a single-panel normalized flux plot (for compact display).
    
    Parameters
    ----------
    rb_spec : rb_spec object
        The spectrum object containing all data and analysis results
    item : object with template and analysis attributes
        Batch item containing metadata
    figure : matplotlib.figure.Figure, optional
        Existing figure to plot into
    clear_figure : bool, optional
        Whether to clear the figure before plotting
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot
    """
    # Create or use existing figure
    if figure is None:
        fig = Figure(figsize=(8, 5))
    else:
        fig = figure
        if clear_figure:
            fig.clear()
    
    # Single subplot
    ax = fig.add_subplot(111)
    
    # Plot normalized flux panel
    _plot_normalized_panel(rb_spec, item, ax, show_xlabel=True)
    
    # Add transition info as title
    transition_name = item.template.transition_name
    transition_wave = item.template.transition
    title = f"{transition_name} ({transition_wave:.2f} Å)"
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # Don't use tight_layout to avoid shrinking
    return fig


def plot_continuum_focus(rb_spec, item, figure=None, clear_figure=True):
    """
    Create a continuum-focused plot showing flux + continuum fit details.
    
    This is useful for continuum editing workflows where you want to focus
    on the continuum fitting quality.
    """
    # Create or use existing figure
    if figure is None:
        fig = Figure(figsize=(10, 6))
    else:
        fig = figure
        if clear_figure:
            fig.clear()
    
    # Single subplot focused on continuum
    ax = fig.add_subplot(111)
    
    # Plot flux, error, and continuum with enhanced detail
    ax.step(rb_spec.velo, rb_spec.flux_slice, 'k-', where='mid', 
           linewidth=0.8, label='Flux', alpha=0.9)
    
    if hasattr(rb_spec, 'error_slice'):
        ax.step(rb_spec.velo, rb_spec.error_slice, 'r-', where='mid', 
               alpha=0.6, linewidth=0.6, label='Error')
    
    # Enhanced continuum display
    if hasattr(rb_spec, 'cont'):
        ax.plot(rb_spec.velo, rb_spec.cont, 'g-', linewidth=2.5, 
               label='Continuum Fit', alpha=0.9)
    
    # Show masked regions more prominently
    if hasattr(item, 'analysis') and hasattr(item.analysis, 'continuum_masks'):
        mask_count = 0
        for vmin, vmax in item.analysis.continuum_masks:
            label = 'Masked Regions' if mask_count == 0 else ""
            ax.axvspan(vmin, vmax, alpha=0.3, color='red', label=label)
            mask_count += 1
    
    # Reference lines
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle=':', alpha=0.7)
    
    # Labels and formatting
    ax.set_xlabel('Velocity (km/s)', fontsize=11)
    ax.set_ylabel('Flux', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Smaller legend
    ax.legend(loc='best', fontsize=8)
    
    # Title with transition info
    transition_name = item.template.transition_name
    title = f"Continuum Fit: {transition_name}"
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # Don't use tight_layout to avoid shrinking
    return fig