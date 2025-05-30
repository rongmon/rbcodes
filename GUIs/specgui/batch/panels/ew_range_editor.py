# rbcodes/GUIs/specgui/batch/panels/ew_range_editor.py
"""
Separate module for EW range editing functionality.
This module provides an enhanced EW range editor with interactive features.
"""

import numpy as np
from datetime import datetime
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                           QSpinBox, QFormLayout, QCheckBox, QDialogButtonBox)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib plots in the EW editor."""
    
    def __init__(self, parent=None, width=4, height=3, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)

def edit_ew_range_dialog(item, current_spec, controller, parent=None):
    """
    Launch an enhanced EW range editor dialog.
    
    Parameters
    ----------
    item : object with template, analysis, results attributes
        The batch item to edit (now includes row_index)
    current_spec : rb_spec
        The current spectrum object
    controller : BatchController
        The batch controller for updates
    parent : QWidget
        Parent widget for the dialog
    
    Returns
    -------
    bool
        True if EW range was updated successfully, False otherwise
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit EW Measurement Range")
    dialog.setMinimumSize(600, 500)
    layout = QVBoxLayout(dialog)
    
    # Get current values
    current_ew_vmin = item.template.ew_vmin
    current_ew_vmax = item.template.ew_vmax
    slice_vmin = item.template.slice_vmin
    slice_vmax = item.template.slice_vmax
    
    # Form for velocity range
    form = QFormLayout()
    
    # Show current slice range for reference
    slice_info = QLabel(f"Slice Range: [{slice_vmin}, {slice_vmax}] km/s")
    slice_info.setStyleSheet("QLabel { color: gray; font-style: italic; }")
    form.addRow("Reference:", slice_info)
    
    # Velocity min/max inputs
    ew_vmin_spin = QSpinBox()
    ew_vmin_spin.setRange(-5000, 0)
    ew_vmin_spin.setValue(current_ew_vmin)
    ew_vmin_spin.setSingleStep(10)
    ew_vmin_spin.setSuffix(" km/s")
    
    ew_vmax_spin = QSpinBox()
    ew_vmax_spin.setRange(0, 5000)
    ew_vmax_spin.setValue(current_ew_vmax)
    ew_vmax_spin.setSingleStep(10)
    ew_vmax_spin.setSuffix(" km/s")
    
    form.addRow("EW Min:", ew_vmin_spin)
    form.addRow("EW Max:", ew_vmax_spin)
    
    # Validation label
    validation_label = QLabel("")
    form.addRow("", validation_label)
    
    # Validation function
    def validate_ew_range():
        ew_min = ew_vmin_spin.value()
        ew_max = ew_vmax_spin.value()
        
        # Check if EW range is valid
        if ew_min >= ew_max:
            validation_label.setText("⚠️ EW min must be less than EW max")
            validation_label.setStyleSheet("QLabel { color: red; }")
            return False
        
        # Check if EW range is within slice range (with tolerance)
        tolerance = 50
        if ew_min < (slice_vmin - tolerance) or ew_max > (slice_vmax + tolerance):
            validation_label.setText(f"⚠️ EW range should be within slice range ± {tolerance} km/s")
            validation_label.setStyleSheet("QLabel { color: red; }")
            return False
        
        validation_label.setText("✓ EW range is valid")
        validation_label.setStyleSheet("QLabel { color: green; }")
        return True
    
    # Connect validation to value changes
    ew_vmin_spin.valueChanged.connect(validate_ew_range)
    ew_vmax_spin.valueChanged.connect(validate_ew_range)
    
    # Initial validation
    validate_ew_range()
    
    # SNR calculation option
    snr_check = QCheckBox("Calculate SNR")
    snr_check.setChecked(item.analysis.calculate_snr)
    
    binsize_spin = QSpinBox()
    binsize_spin.setRange(1, 10)
    binsize_spin.setValue(item.analysis.binsize)
    binsize_spin.setEnabled(snr_check.isChecked())
    
    snr_check.toggled.connect(binsize_spin.setEnabled)
    
    snr_layout = QHBoxLayout()
    snr_layout.addWidget(snr_check)
    snr_layout.addWidget(QLabel("Bin Size:"))
    snr_layout.addWidget(binsize_spin)
    
    form.addRow("SNR Options:", snr_layout)
    
    # Interactive selection option
    interactive_selection = QCheckBox("Enable Interactive Selection")
    interactive_selection.setChecked(False)
    interactive_selection.setToolTip("Enable clicking on the plot to set velocity limits")
    form.addRow("", interactive_selection)
    
    # Tooltip help
    help_text = QLabel("When enabled, click on plot to set integration limits. Press 'r' to reset.")
    help_text.setWordWrap(True)
    help_text.setStyleSheet("color: gray; font-style: italic;")
    form.addRow("", help_text)
    
    layout.addLayout(form)
    
    # Add preview plot if spectrum is available
    if current_spec and hasattr(current_spec, 'velo'):
        preview_label = QLabel("Current measurement range:")
        layout.addWidget(preview_label)
        
        canvas = MatplotlibCanvas(dialog, width=5, height=3, dpi=100)
        spec_obj = current_spec
        
        # Initial plot
        canvas.axes.step(spec_obj.velo, spec_obj.fnorm, 'k-', where='mid')
        canvas.axes.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
        canvas.axes.axvspan(current_ew_vmin, current_ew_vmax, alpha=0.2, color='blue')
        canvas.axes.set_xlabel('Velocity (km/s)')
        canvas.axes.set_ylabel('Normalized Flux')
        canvas.axes.set_title('Current EW Range')
        canvas.draw()
        
        # Function to update the preview plot
        def update_preview():
            canvas.axes.clear()
            canvas.axes.step(spec_obj.velo, spec_obj.fnorm, 'k-', where='mid')
            canvas.axes.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
            canvas.axes.axvspan(ew_vmin_spin.value(), ew_vmax_spin.value(), alpha=0.2, color='blue')
            canvas.axes.set_xlabel('Velocity (km/s)')
            canvas.axes.set_ylabel('Normalized Flux')
            canvas.axes.set_title('Updated EW Range')
            
            # Add markers at vmin and vmax if interactive selection is enabled
            if interactive_selection.isChecked():
                vmin_val = ew_vmin_spin.value()
                vmax_val = ew_vmax_spin.value()
                canvas.axes.plot([vmin_val], [1.0], 'bo', ms=6, zorder=10)
                canvas.axes.plot([vmax_val], [1.0], 'bo', ms=6, zorder=10)
                
                canvas.axes.text(vmin_val, 1.05, f'vmin: {vmin_val}', 
                                ha='center', va='bottom', fontsize=8, color='blue')
                canvas.axes.text(vmax_val, 1.05, f'vmax: {vmax_val}', 
                                ha='center', va='bottom', fontsize=8, color='blue')
            
            canvas.draw()
        
        # Interactive selection event handlers
        def on_canvas_click(event):
            if not event.inaxes or not interactive_selection.isChecked():
                return
                
            velocity = event.xdata
            current_vmin = ew_vmin_spin.value()
            current_vmax = ew_vmax_spin.value()
            
            if velocity < current_vmin:
                ew_vmin_spin.setValue(int(velocity))
            elif velocity > current_vmax:
                ew_vmax_spin.setValue(int(velocity))
            else:
                if abs(velocity - current_vmin) < abs(velocity - current_vmax):
                    ew_vmin_spin.setValue(int(velocity))
                else:
                    ew_vmax_spin.setValue(int(velocity))
        
        def on_canvas_key_press(event):
            if not interactive_selection.isChecked():
                return
            if event.key == 'r':
                ew_vmin_spin.setValue(current_ew_vmin)
                ew_vmax_spin.setValue(current_ew_vmax)
        
        def toggle_interactive_mode(enabled):
            if enabled:
                canvas.mpl_connect('button_press_event', on_canvas_click)
                canvas.mpl_connect('key_press_event', on_canvas_key_press)
                canvas.setFocus()
                update_preview()
            else:
                update_preview()
        
        interactive_selection.toggled.connect(toggle_interactive_mode)
        ew_vmin_spin.valueChanged.connect(update_preview)
        ew_vmax_spin.valueChanged.connect(update_preview)
        
        layout.addWidget(canvas)
    
    # Add buttons
    button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    
    def accept_with_validation():
        if not validate_ew_range():
            controller.status_updated.emit("Please fix the EW range validation errors before proceeding.")
            return
        dialog.accept()
    
    button_box.accepted.connect(accept_with_validation)
    button_box.rejected.connect(dialog.reject)
    layout.addWidget(button_box)
    
    # Show dialog
    if dialog.exec_() == QDialog.Accepted:
        # Get new values
        new_ew_vmin = ew_vmin_spin.value()
        new_ew_vmax = ew_vmax_spin.value()
        new_calculate_snr = snr_check.isChecked()
        new_binsize = binsize_spin.value()
        
        # Update EW range using row index
        success = _update_ew_range_in_master_table(
            item, 
            new_ew_vmin, 
            new_ew_vmax, 
            new_calculate_snr, 
            new_binsize,
            controller
        )
        
        return success
    
    return False

def _update_ew_range_in_master_table(item, new_ew_vmin, new_ew_vmax, calculate_snr, binsize, controller):
    """Update EW range and recalculate measurements, ensuring master table consistency."""
    try:
        # Get row index from item
        row_index = getattr(item, 'row_index', None)
        if row_index is None:
            print("ERROR: No row_index found in item")
            return False
        
        
        # Update template parameters in master table
        controller.master_table.update_template(row_index, 
                                               ew_vmin=new_ew_vmin, 
                                               ew_vmax=new_ew_vmax)
        
        # Update analysis settings in master table
        controller.master_table.update_analysis(row_index,
                                               calculate_snr=calculate_snr,
                                               binsize=binsize,
                                               last_modified=datetime.now().isoformat())
        
        # Generate fresh spectrum with new parameters
        updated_item = controller.master_table.get_item(row_index)
        if not updated_item:
            print("ERROR: Failed to get updated item from master table")
            return False
        
        updated_spec = _generate_spectrum_from_item(updated_item)
        if not updated_spec:
            print("ERROR: Failed to generate spectrum")
            return False
        
        # Recalculate EW with new range
        updated_spec.compute_EW(
            updated_item.template.transition,
            vmin=new_ew_vmin,
            vmax=new_ew_vmax,
            SNR=calculate_snr,
            _binsize=binsize
        )
        
        
        # Update results in master table
        controller.master_table.update_results(
            row_index,
            W=updated_spec.W,
            W_e=updated_spec.W_e,
            N=updated_spec.N,
            N_e=updated_spec.N_e,
            logN=updated_spec.logN,
            logN_e=updated_spec.logN_e,
            vel_centroid=getattr(updated_spec, 'vel_centroid', 0),
            vel_disp=getattr(updated_spec, 'vel_disp', 0),
            SNR=getattr(updated_spec, 'SNR', 0),
            calculation_timestamp=datetime.now().isoformat()
        )
        
        # Mark as complete
        controller.master_table.update_analysis(row_index, processing_status="complete")
        
        # Store the updated rb_spec object in master table
        controller.master_table.set_rb_spec_object(row_index, updated_spec)
        
        return True
        
    except Exception as e:
        print(f"ERROR: EW range update failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def _generate_spectrum_from_item(item):
    """Generate rb_spec object from item parameters (local copy of the function)."""
    try:
        from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object
        import numpy as np
        import os
        
        # Load spectrum file
        if item.template.filename.lower().endswith('.json'):
            spec = load_rb_spec_object(item.template.filename)
        else:
            # Determine file type
            ext = os.path.splitext(item.template.filename)[1].lower()
            filetype = 'linetools' if ext == '.fits' else ('ascii' if ext in ['.txt', '.dat'] else None)
            spec = rb_spec.from_file(item.template.filename, filetype=filetype)
        
        # Apply saved parameters
        spec.shift_spec(item.template.redshift)
        spec.slice_spec(
            item.template.transition,
            item.template.slice_vmin, item.template.slice_vmax,
            use_vel=True,
            linelist=item.template.linelist
        )
        
        # Apply continuum (from saved parameters)
        if item.analysis.continuum_method == 'polynomial':
            mask = []
            for vmin, vmax in item.analysis.continuum_masks:
                mask.extend([vmin, vmax])
            spec.fit_continuum(
                mask=mask if mask else False,
                Legendre=item.analysis.continuum_order,
                use_weights=item.analysis.use_weights
            )
        else:  # flat
            spec.cont = np.ones_like(spec.flux_slice)
            spec.fnorm = spec.flux_slice.copy()
            spec.enorm = spec.error_slice.copy()
        
        # Apply EW calculation results if we have them
        if item.results.W > 0:
            spec.W = item.results.W
            spec.W_e = item.results.W_e
            spec.N = item.results.N
            spec.N_e = item.results.N_e
            spec.logN = item.results.logN
            spec.logN_e = item.results.logN_e
            spec.vel_centroid = item.results.vel_centroid
            spec.vel_disp = item.results.vel_disp
            spec.SNR = item.results.SNR
            spec.vmin = item.template.ew_vmin
            spec.vmax = item.template.ew_vmax
            spec.trans = item.template.transition_name
            spec.trans_wave = item.template.transition
        
        return spec
        
    except Exception as e:
        print(f"Error generating spectrum: {e}")
        import traceback
        traceback.print_exc()
        return None