# rbcodes/GUIs/specgui/batch/batch_figure_generator.py
"""
Batch figure generation module for rb_spec batch processing.

This module provides functionality to create publication-quality figures
from batch processing results, supporting both individual figures and
multi-page PDF layouts.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from datetime import datetime

# Set matplotlib parameters for consistent styling
plt.rcParams.update({
    'font.size': 8,
    'axes.labelsize': 8,
    'axes.titlesize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 200,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5
})

def create_individual_figures(items, output_directory, file_format='pdf', master_table=None):
    """
    Create individual figure files for each batch item.
    
    Parameters
    ----------
    items : list
        List of batch items from master table
    output_directory : str
        Directory to save figures
    file_format : str
        Output format ('pdf' or 'png')
    master_table : MasterBatchTable
        Master table containing rb_spec objects
    
    Returns
    -------
    success_count : int
        Number of figures successfully created
    error_count : int
        Number of figures that failed to create
    """
    success_count = 0
    error_count = 0
    
    for i, item in enumerate(items):
        try:
            # Get rb_spec object using row index
            spec_object = _get_spec_object(item, master_table, i)
            if spec_object is None:
                print(f"Cannot create figure for {item.template.transition_name}: No rb_spec object")
                error_count += 1
                continue
            
            # Generate filename
            if item.template.filename.lower().endswith('.json'):
                # For JSON files, use the exact same basename
                figure_filename = f"{_get_clean_basename(item)}.{file_format}"
            else:
                # For other files, use the descriptive naming convention
                basename = _get_clean_basename(item)
                transition_id = f"{item.template.transition_name}_{item.template.transition:.0f}"
                figure_filename = f"{basename}_{transition_id}_z{item.template.redshift:.3f}.{file_format}"
            figure_path = os.path.join(output_directory, figure_filename)
            
            # Create individual figure
            fig = _create_single_figure(item, spec_object)
            
            # Save figure
            fig.savefig(figure_path, dpi=200, bbox_inches='tight')
            plt.close(fig)
            
            print(f"Created figure: {figure_filename}")
            success_count += 1
            
        except Exception as e:
            print(f"Error creating figure for {item.template.transition_name}: {str(e)}")
            error_count += 1
    
    return success_count, error_count

def create_multipage_pdf(items, output_path, master_table=None):
    """
    Create a multi-page PDF with 6 transitions per page.
    
    Parameters
    ----------
    items : list
        List of batch items from master table
    output_path : str
        Path for output PDF file
    master_table : MasterBatchTable
        Master table containing rb_spec objects
    
    Returns
    -------
    success : bool
        Whether the PDF was created successfully
    """
    try:
        # Filter items that have rb_spec objects
        valid_items = []
        for i, item in enumerate(items):
            spec_object = _get_spec_object(item, master_table, i)
            if spec_object is not None:
                valid_items.append((item, spec_object))
        
        if not valid_items:
            print("No valid items with rb_spec objects found")
            return False
        
        # Create multi-page PDF
        with PdfPages(output_path) as pdf:
            # Process items in chunks of 6
            for i in range(0, len(valid_items), 6):
                chunk = valid_items[i:i+6]
                
                # Create page with up to 6 transitions
                fig = _create_multipage_figure(chunk, i//6 + 1, len(valid_items))
                
                # Save page to PDF
                pdf.savefig(fig, dpi=200, bbox_inches='tight')
                plt.close(fig)
        
        print(f"Created multi-page PDF: {os.path.basename(output_path)}")
        print(f"Total items: {len(valid_items)}, Total pages: {(len(valid_items) + 5) // 6}")
        return True
        
    except Exception as e:
        print(f"Error creating multi-page PDF: {str(e)}")
        return False

def _get_spec_object(item, master_table, row_index):
    """Get rb_spec object - use existing one from master table."""
    
    if master_table is None:
        print(f"No master table provided for {item.template.transition_name}")
        return None
    
    # Get existing object from master table
    spec_object = master_table.get_rb_spec_object(row_index)
    
    if spec_object is not None:
        return spec_object
    
    # If not found, object wasn't processed yet
    print(f"No rb_spec object found for {item.template.transition_name}")
    return None

def export_batch_figures(items, output_directory, file_format='pdf', figure_type='individual', master_table=None):
    """
    Main function to export batch figures.
    
    Parameters
    ----------
    items : list
        List of batch items from master table
    output_directory : str
        Directory to save figures
    file_format : str
        Output format ('pdf' or 'png')
    figure_type : str
        Type of output ('individual' or 'multipage')
    master_table : MasterBatchTable
        Master table containing rb_spec objects
    
    Returns
    -------
    success : bool
        Whether the export was successful
    message : str
        Status message
    """
    try:
        # Create output directory if it doesn't exist
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        
        if figure_type == 'individual':
            success_count, error_count = create_individual_figures(items, output_directory, file_format, master_table)
            
            if error_count == 0:
                return True, f"Successfully created {success_count} individual figures"
            else:
                return False, f"Created {success_count} figures, {error_count} failed"
        
        elif figure_type == 'multipage' and file_format == 'pdf':
            # Generate filename for multi-page PDF
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            pdf_filename = f"batch_results_{timestamp}.pdf"
            pdf_path = os.path.join(output_directory, pdf_filename)
            
            success = create_multipage_pdf(items, pdf_path, master_table)
            
            if success:
                return True, f"Successfully created multi-page PDF: {pdf_filename}"
            else:
                return False, "Failed to create multi-page PDF"
        
        else:
            return False, "Multi-page option is only available for PDF format"
    
    except Exception as e:
        return False, f"Error during figure export: {str(e)}"


def _get_clean_basename(item):
    """Get clean basename for filename generation."""
    # Simple logic: if JSON input, use the JSON name; otherwise use file name
    if item.template.filename.lower().endswith('.json'):
        # For JSON files, use the exact same basename (no extension)
        return os.path.splitext(os.path.basename(item.template.filename))[0]
    else:
        # For other files, use the standard approach
        return os.path.splitext(os.path.basename(item.template.filename))[0]

def _create_single_figure(item, spec_object):
    """Create a two-panel figure for a single transition."""
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    
    # Left panel: Flux + continuum
    _plot_flux_panel(spec_object, item, ax1)
    
    # Right panel: Normalized flux + EW region
    _plot_normalized_panel(spec_object, item, ax2)
    
    # Add overall title
    title = f"{item.template.transition_name} ({item.template.transition:.2f} Å) at z = {item.template.redshift:.6f}"
    fig.suptitle(title, fontsize=10, y=0.95)
    
    plt.tight_layout()
    return fig

def _create_multipage_figure(chunk, page_num, total_items):
    """Create a figure with up to 6 transitions (6 rows × 2 columns layout)."""
    n_items = len(chunk)
    
    # Create figure with 6 rows, 2 columns
    fig, axes = plt.subplots(6, 2, figsize=(12, 16))
    
    # Add page title
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    page_title = f"Batch Processing Results - {timestamp}"
    fig.suptitle(page_title, fontsize=12, y=0.98)
    
    # Add page number
    total_pages = (total_items + 5) // 6
    fig.text(0.5, 0.02, f"Page {page_num} of {total_pages}", ha='center', fontsize=8)
    
    # Plot each transition
    for i, (item, spec_object) in enumerate(chunk):
        row = i  # Items go in rows 0, 1, 2, 3, 4, 5
        
        # Left column: flux panel
        ax_flux = axes[row, 0]
        _plot_flux_panel(spec_object, item, ax_flux, show_xlabel=(row == 5 or row == n_items-1))
        
        # Add transition title inside the left panel
        title = f"{item.template.transition_name}, z = {item.template.redshift:.3f}"
        ax_flux.text(0.02, 0.95, title, transform=ax_flux.transAxes, 
                    fontsize=9, weight='bold', verticalalignment='top',
                    bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
        
        # Right column: normalized panel  
        ax_norm = axes[row, 1]
        _plot_normalized_panel(spec_object, item, ax_norm, show_xlabel=(row == 5 or row == n_items-1))
    
    # Hide unused rows
    for i in range(n_items, 6):
        axes[i, 0].set_visible(False)
        axes[i, 1].set_visible(False)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.08, right=0.98, hspace=0.3)
    
    return fig

def _plot_flux_panel(spec_object, item, ax, show_xlabel=True):
    """Plot flux panel with continuum and masked regions."""
    # Plot flux and error
    ax.step(spec_object.velo, spec_object.flux_slice, 'k-', linewidth=0.5, where='mid', label='Flux')
    ax.step(spec_object.velo, spec_object.error_slice, 'r-', linewidth=0.5, where='mid', alpha=0.7, label='Error')
    
    # Plot continuum
    ax.plot(spec_object.velo, spec_object.cont, color='purple', linewidth=1.0, label='Continuum')
    
    # Add masked regions
    for vmin, vmax in item.analysis.continuum_masks:
        ax.axvspan(vmin, vmax, alpha=0.2, color='red', label='Masked' if vmin == item.analysis.continuum_masks[0][0] else "")
    
    # Reference lines
    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(x=0, color='black', linestyle=':', linewidth=0.5, alpha=0.5)
    
    # Labels and formatting
    ax.set_ylabel('Flux')
    if show_xlabel:
        ax.set_xlabel('Velocity (km/s)')
    
    # Set reasonable limits
    velocity = spec_object.velo
    flux = spec_object.flux_slice
    
    x_buffer = (max(velocity) - min(velocity)) * 0.05
    ax.set_xlim(min(velocity) - x_buffer, max(velocity) + x_buffer)
    
    y_min = min(0, min(flux)) - 0.1 * abs(max(flux))
    y_max = max(flux) + 0.2 * abs(max(flux))
    ax.set_ylim(y_min, y_max)
    
    ax.grid(True, alpha=0.3)

def _plot_normalized_panel(spec_object, item, ax, show_xlabel=True):
    """Plot normalized flux panel with EW region and results."""
    # Plot normalized flux and error
    ax.step(spec_object.velo, spec_object.fnorm, 'k-', linewidth=0.5, where='mid', label='Normalized Flux')
    ax.step(spec_object.velo, spec_object.enorm, 'r-', linewidth=0.5, where='mid', alpha=0.7, label='Error')
    
    # Add EW integration region
    ew_vmin = item.template.ew_vmin
    ew_vmax = item.template.ew_vmax
    ax.axvspan(ew_vmin, ew_vmax, alpha=0.15, color='blue', label='EW Region')
    ax.axvline(x=ew_vmin, color='blue', linestyle=':', linewidth=0.5)
    ax.axvline(x=ew_vmax, color='blue', linestyle=':', linewidth=0.5)
    
    # Reference lines
    ax.axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axhline(y=1, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(x=0, color='black', linestyle=':', linewidth=0.5, alpha=0.5)
    
    # Add results text box
    ew = item.results.W
    ew_err = item.results.W_e
    logn = item.results.logN
    logn_err = item.results.logN_e
    snr = item.results.SNR

    # Handle negative N case for figure display
    if item.results.N < 0 and item.results.N_e > 0:
        import numpy as np
        logn = 0.0
        logn_err = np.log10(item.results.N_e)
    
    results_text = f"EW = {ew:.3f} ± {ew_err:.3f} Å\nlog N = {logn:.2f} ± {logn_err:.2f}\nSNR = {snr:.1f}"
    
    ax.text(0.98, 0.05, results_text, transform=ax.transAxes, 
            horizontalalignment='right', verticalalignment='bottom',
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'),
            fontsize=7)
    
    # Labels and formatting
    ax.set_ylabel('Normalized Flux')
    if show_xlabel:
        ax.set_xlabel('Velocity (km/s)')
    
    # Set reasonable limits
    velocity = spec_object.velo
    norm_flux = spec_object.fnorm
    
    x_buffer = (max(velocity) - min(velocity)) * 0.05
    ax.set_xlim(min(velocity) - x_buffer, max(velocity) + x_buffer)
    
    # Set y-limits for normalized flux
    y_min = max(-0.1, min(norm_flux) - 0.1)
    y_max = min(2.0, max(norm_flux) + 0.3)
    ax.set_ylim(y_min, y_max)
    
    ax.grid(True, alpha=0.3)

