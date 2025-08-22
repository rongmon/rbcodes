import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                           QLabel, QLineEdit, QPushButton, QFileDialog, QTabWidget, 
                           QGridLayout, QGroupBox, QFormLayout, QDoubleSpinBox, QSpinBox, 
                           QCheckBox, QMessageBox, QTableWidget, QTableWidgetItem, QHeaderView,
                           QSplitter, QFrame, QComboBox, QSizePolicy, QMenu, QAction,
                           QProgressBar, QStatusBar, QTextEdit, QDialog,QScrollArea)
from PyQt5.QtCore import Qt, QSettings, QTimer
from PyQt5.QtGui import QIcon, QFont, QTextCursor
import datetime  # Add this to fix the datetime error

# Import the LLSFitter class
from rbcodes.IGM.LLSFitter import LLSFitter

class MatplotlibWidget(QWidget):
    def __init__(self, parent=None):
        super(MatplotlibWidget, self).__init__(parent)
        
        # Create figure and canvas
        self.figure = plt.figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        # Set up the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
    
    def clear_figure(self):
        self.figure.clear()
        self.canvas.draw()
    
    def get_figure(self):
        return self.figure


class ContinuumRegionsTable(QTableWidget):
    def __init__(self, parent=None):
        super(ContinuumRegionsTable, self).__init__(parent)
        
        # Set up the table
        self.setColumnCount(2)
        self.setHorizontalHeaderLabels(["Min λ (Å)", "Max λ (Å)"])
        
        # Set column widths
        header = self.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        
        # Add default regions
        #self.default_regions = [
        #    (860, 872), (893, 899), (897, 899), (905, 910),
        #    (940, 944), (946, 948), (927, 928), (918.5, 919),
        #    (931, 933), (934, 936), (951, 970)
        #]

        self.default_regions = [
             (893, 899), (897, 899), (905, 910),
            (940, 944), (946, 948), (927, 928), (918.5, 919),
            (931, 933), (934, 936), (951, 970)
        ]

        
        # Connect signals
        self.cellChanged.connect(self.validate_cell)
    
    def set_default_regions(self):
        """Set the table to contain the default regions"""
        self.blockSignals(True)  # Block signals to prevent recursive calls
        self.setRowCount(len(self.default_regions))
        
        for i, (wmin, wmax) in enumerate(self.default_regions):
            self.setItem(i, 0, QTableWidgetItem(str(wmin)))
            self.setItem(i, 1, QTableWidgetItem(str(wmax)))
        
        self.blockSignals(False)  # Unblock signals
    
    def get_regions(self):
        """Get the current regions as a list of tuples"""
        regions = []
        for i in range(self.rowCount()):
            try:
                wmin = float(self.item(i, 0).text())
                wmax = float(self.item(i, 1).text())
                regions.append((wmin, wmax))
            except (ValueError, AttributeError):
                # Skip invalid rows
                pass
        
        return regions
    
    def validate_cell(self, row, column):
        """Validate that cell contents are numeric"""
        try:
            item = self.item(row, column)
            if item is not None:
                value = float(item.text())
                
                # Check if min < max for the row
                if column == 0 and self.item(row, 1) is not None and self.item(row, 1) != 0:
                    wmax = float(self.item(row, 1).text())
                    # Modified to ignore when max is 0 (unset value)
                    if value >= wmax and wmax != 0.0:
                        QMessageBox.warning(self, "Invalid Value", 
                                          "Minimum wavelength must be less than maximum wavelength")
                        item.setText("0.0")
                
                if column == 1 and self.item(row, 0) is not None and self.item(row, 0) != 0:
                    wmin = float(self.item(row, 0).text())
                    # Modified to ignore when min is 0 (unset value)
                    if value <= wmin and wmin != 0.0:
                        QMessageBox.warning(self, "Invalid Value", 
                                          "Maximum wavelength must be greater than minimum wavelength")
                        item.setText("0.0")
            
        except ValueError:
            QMessageBox.warning(self, "Invalid Value", "Please enter a valid number")
            item.setText("0.0")

    
    # For the ContinuumRegionsTable class, add a sort_regions method and modify add_row to call it
    
    def sort_regions(self):
        """Sort the continuum regions by minimum wavelength"""
        # Get all regions
        regions = self.get_regions()
        if not regions:
            return
        
        # Sort regions by minimum wavelength
        regions.sort(key=lambda r: r[0])
        
        # Clear the table
        self.blockSignals(True)  # Block signals to prevent recursive calls
        self.setRowCount(0)
        
        # Add sorted regions back
        for wmin, wmax in regions:
            row = self.rowCount()
            self.insertRow(row)
            self.setItem(row, 0, QTableWidgetItem(str(wmin)))
            self.setItem(row, 1, QTableWidgetItem(str(wmax)))
        
        self.blockSignals(False)  # Unblock signals
    
    def add_row(self):
        """Add a new empty row to the table and sort"""
        row = self.rowCount()
        self.insertRow(row)
        self.setItem(row, 0, QTableWidgetItem("0.0"))
        self.setItem(row, 1, QTableWidgetItem("0.0"))
        self.sort_regions()  # Sort after adding    
    
    def remove_selected_rows(self):
        """Remove the selected rows"""
        rows = set()
        for item in self.selectedItems():
            rows.add(item.row())
        
        # Remove in reverse order to avoid index shift
        for row in sorted(rows, reverse=True):
            self.removeRow(row)


class ResultsDialog(QDialog):
    """Dialog to display detailed fit results"""
    def __init__(self, results, parent=None):
        super(ResultsDialog, self).__init__(parent)
        self.setWindowTitle("Detailed Fit Results")
        self.resize(600, 400)
        
        layout = QVBoxLayout(self)
        
        # Create text display
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        font = QFont("Monospace")
        font.setStyleHint(QFont.TypeWriter)
        text_edit.setFont(font)
        
        # Format and display results
        text = "LLS Fit Results\n"
        text += "==============\n\n"
        
        if 'curve_fit' in results:
            cf = results['curve_fit']
            text += "Curve Fit Results:\n"
            text += "-----------------\n"
            text += f"C0:      {cf['C0']:.6f} ± {cf['C0_err']:.6f}\n"
            text += f"C1:      {cf['C1']:.6f} ± {cf['C1_err']:.6f}\n"
            text += f"log N(HI): {cf['logNHI']:.4f} ± {cf['logNHI_err']:.4f}\n\n"
        
        if 'mcmc' in results:
            mc = results['mcmc']
            text += "MCMC Results:\n"
            text += "------------\n"
            text += f"C0:      {mc['C0']:.6f} ± {mc['C0_err']:.6f}\n"
            text += f"C1:      {mc['C1']:.6f} ± {mc['C1_err']:.6f}\n"
            text += f"log N(HI): {mc['logNHI']:.4f} ± {mc['logNHI_err']:.4f}\n\n"
            
            # Add percentiles if available
            if 'samples' in results:
                samples = results['samples']
                if samples is not None:
                    text += "MCMC Percentiles:\n"
                    text += "----------------\n"
                    for i, param in enumerate(['C0', 'C1', 'log N(HI)']):
                        p16, p50, p84 = np.percentile(samples[:, i], [16, 50, 84])
                        text += f"{param}: {p50:.6f} (+{p84-p50:.6f}) (-{p50-p16:.6f})\n"
        
        text_edit.setText(text)
        layout.addWidget(text_edit)
        
        # Add buttons
        button_layout = QHBoxLayout()
        copy_button = QPushButton("Copy to Clipboard")
        copy_button.clicked.connect(lambda: self.copy_to_clipboard(text))
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        
        button_layout.addWidget(copy_button)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
    
    def copy_to_clipboard(self, text):
        """Copy results to clipboard"""
        clipboard = QApplication.clipboard()
        clipboard.setText(text)
        QMessageBox.information(self, "Copied", "Results copied to clipboard")


class LLSFitterGUI(QMainWindow):
    def __init__(self):
        super(LLSFitterGUI, self).__init__()
        
        # Initialize LLSFitter object to None (will create when loading a spectrum)
        self.lls_fitter = None
        
        # Set up the UI
        self.init_ui()
        
        # Create menus
        self.create_menus()
        
        # Load settings
        self.settings = QSettings("AstronomyTools", "LLSFitterGUI")
        self.load_settings()
        
        # Set window properties
        self.setWindowTitle("LLS Fitter GUI")
        self.resize(1200, 800)
        
    def save_continuum_regions(self):
        """Save current continuum regions to a template file"""
        regions = self.get_current_continuum_regions()
        if not regions:
            QMessageBox.warning(self, "No Regions", "No valid continuum regions to save")
            return
        
        # Open file dialog
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Continuum Regions", "", "JSON Files (*.json);;All Files (*.*)"
        )
        
        if not file_path:
            return
        
        try:
            import json
            
            # Create data structure with metadata
            data = {
                "regions": regions,
                "metadata": {
                    "created": datetime.datetime.now().isoformat(),
                    "description": "LLS Fitter continuum regions template"
                }
            }
            
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2)
            
            self.status_bar.showMessage(f"Continuum regions saved to {file_path}", 3000)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error saving continuum regions: {str(e)}")
    
    def load_continuum_regions(self):
        """Load continuum regions from a template file"""
        # Open file dialog
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Continuum Regions", "", "JSON Files (*.json);;All Files (*.*)"
        )
        
        if not file_path:
            return
        
        try:
            import json
            
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            if "regions" not in data:
                QMessageBox.warning(self, "Invalid File", "The file does not contain valid continuum regions")
                return
            
            regions = data["regions"]
            
            # Clear current table
            self.continuum_table.setRowCount(0)
            
            # Add loaded regions
            for wmin, wmax in regions:
                row = self.continuum_table.rowCount()
                self.continuum_table.insertRow(row)
                self.continuum_table.setItem(row, 0, QTableWidgetItem(str(wmin)))
                self.continuum_table.setItem(row, 1, QTableWidgetItem(str(wmax)))
            
            self.status_bar.showMessage(f"Loaded {len(regions)} continuum regions from {file_path}", 3000)
            
            # Preview the loaded regions
            self.preview_continuum_regions()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading continuum regions: {str(e)}")
    
    def suggest_continuum_regions(self):
        """Suggest continuum regions based on the spectrum"""
        if self.lls_fitter is None:
            QMessageBox.warning(self, "Error", "No spectrum loaded")
            return
        
        try:
            # A simple algorithm to suggest continuum regions
            # This is just a basic implementation - in practice, this would be more sophisticated
            
            # Get the spectrum data
            wv = self.lls_fitter.rest_wave
            flux = self.lls_fitter.flux
            
            # Define wavelength range to analyze
            wmin = 850
            wmax = 1000
            
            # Create mask for the wavelength range
            mask = (wv >= wmin) & (wv <= wmax)
            
            # Skip if no points in range
            if not np.any(mask):
                QMessageBox.warning(self, "Error", "No data points in the analysis range")
                return
            
            # Calculate median and standard deviation of flux
            median_flux = np.median(flux[mask])
            std_flux = np.std(flux[mask])
            
            # Set threshold for continuum (this is a simple approach)
            lower_threshold = median_flux - 0.5 * std_flux
            upper_threshold = median_flux + 2.0 * std_flux
            
            # Find regions where flux is within thresholds
            continuum_mask = (flux >= lower_threshold) & (flux <= upper_threshold) & mask
            
            # Group adjacent points into regions
            regions = []
            in_region = False
            region_start = None
            
            for i in range(len(wv)):
                if continuum_mask[i] and not in_region:
                    # Start of new region
                    in_region = True
                    region_start = wv[i]
                elif not continuum_mask[i] and in_region:
                    # End of region
                    in_region = False
                    region_end = wv[i-1]
                    
                    # Only add if region is wide enough
                    if region_end - region_start >= 2.0:  # Minimum width in Angstroms
                        regions.append((float(region_start), float(region_end)))
            
            # Handle case where last point is part of a region
            if in_region:
                region_end = wv[-1]
                if region_end - region_start >= 2.0:
                    regions.append((float(region_start), float(region_end)))
            
            # Limit to a reasonable number of regions
            max_regions = 15
            if len(regions) > max_regions:
                # Sort by width and keep the widest
                regions.sort(key=lambda r: r[1] - r[0], reverse=True)
                regions = regions[:max_regions]
            
            # Update table with suggested regions
            if regions:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Question)
                msg.setText(f"Found {len(regions)} potential continuum regions. Would you like to:")
                msg.setWindowTitle("Suggested Continuum Regions")
                replace_button = msg.addButton("Replace Current", QMessageBox.ActionRole)
                add_button = msg.addButton("Add to Current", QMessageBox.ActionRole)
                cancel_button = msg.addButton("Cancel", QMessageBox.RejectRole)
                
                msg.exec_()
                
                if msg.clickedButton() == replace_button:
                    # Clear current table
                    self.continuum_table.setRowCount(0)
                    
                    # Add suggested regions
                    for wmin, wmax in regions:
                        row = self.continuum_table.rowCount()
                        self.continuum_table.insertRow(row)
                        self.continuum_table.setItem(row, 0, QTableWidgetItem(str(wmin)))
                        self.continuum_table.setItem(row, 1, QTableWidgetItem(str(wmax)))
                    
                    self.status_bar.showMessage(f"Added {len(regions)} suggested continuum regions", 3000)
                    
                elif msg.clickedButton() == add_button:
                    # Add suggested regions to current ones
                    for wmin, wmax in regions:
                        row = self.continuum_table.rowCount()
                        self.continuum_table.insertRow(row)
                        self.continuum_table.setItem(row, 0, QTableWidgetItem(str(wmin)))
                        self.continuum_table.setItem(row, 1, QTableWidgetItem(str(wmax)))
                    
                    self.status_bar.showMessage(f"Added {len(regions)} suggested continuum regions", 3000)
                
                # Preview if regions were added
                if msg.clickedButton() != cancel_button:
                    self.preview_continuum_regions()
            else:
                QMessageBox.information(self, "No Regions Found", 
                                      "Could not identify suitable continuum regions automatically")
        
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error suggesting continuum regions: {str(e)}")
            
    def create_menus(self):
        """Create application menus"""
        # Create menu bar
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("&File")
        
        open_action = QAction("&Open Spectrum...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.browse_file)
        file_menu.addAction(open_action)
        
        file_menu.addSeparator()
        
        save_results_action = QAction("&Save Results...", self)
        save_results_action.setShortcut("Ctrl+S")
        save_results_action.triggered.connect(self.save_results)
        file_menu.addAction(save_results_action)

        load_results_action = QAction("&Load Results...", self)
        load_results_action.setShortcut("Ctrl+L")
        load_results_action.triggered.connect(self.load_results)
        file_menu.addAction(load_results_action)
        
        export_plot_action = QAction("&Export Current Plot...", self)
        export_plot_action.setShortcut("Ctrl+E")
        export_plot_action.triggered.connect(self.export_current_plot)
        file_menu.addAction(export_plot_action)
        
        file_menu.addSeparator()
        
        # Continuum regions submenu
        continuum_menu = file_menu.addMenu("Continuum &Regions")
        
        save_regions_action = QAction("&Save Regions Template...", self)
        save_regions_action.triggered.connect(self.save_continuum_regions)
        continuum_menu.addAction(save_regions_action)
        
        load_regions_action = QAction("&Load Regions Template...", self)
        load_regions_action.triggered.connect(self.load_continuum_regions)
        continuum_menu.addAction(load_regions_action)
        
        suggest_regions_action = QAction("&Suggest Regions...", self)
        suggest_regions_action.triggered.connect(self.suggest_continuum_regions)
        continuum_menu.addAction(suggest_regions_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction("E&xit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Analysis menu
        analysis_menu = menubar.addMenu("&Analysis")
        
        preview_action = QAction("&Preview Continuum Regions", self)
        preview_action.triggered.connect(self.preview_continuum_regions)
        analysis_menu.addAction(preview_action)
        
        analysis_menu.addSeparator()
        
        curve_fit_action = QAction("Run &Curve Fit", self)
        curve_fit_action.triggered.connect(self.run_curve_fit)
        analysis_menu.addAction(curve_fit_action)
        
        mcmc_fit_action = QAction("Run &MCMC Fit", self)
        mcmc_fit_action.triggered.connect(self.run_mcmc_fit)
        analysis_menu.addAction(mcmc_fit_action)
        
        analysis_menu.addSeparator()
        
        view_results_action = QAction("&View Detailed Results", self)
        view_results_action.triggered.connect(self.view_detailed_results)
        analysis_menu.addAction(view_results_action)

        save_posterior_action = QAction("Save MCMC Posterior...", self)
        save_posterior_action.triggered.connect(self.save_mcmc_posterior)
        analysis_menu.addAction(save_posterior_action)
        
        # Help menu
        help_menu = menubar.addMenu("&Help")
        
        about_action = QAction("&About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
        
        help_action = QAction("&Help", self)
        help_action.triggered.connect(self.show_help)
        help_menu.addAction(help_action)
    
    def save_mcmc_posterior(self):
        """Save the MCMC posterior samples to a file"""
        if self.lls_fitter is None or self.lls_fitter.mcmc_results is None:
            QMessageBox.warning(self, "No Results", "No MCMC results available to save")
            return
        
        # Open file dialog for CSV or NPY
        file_path, selected_filter = QFileDialog.getSaveFileName(
            self, "Save MCMC Posterior", "", 
            "CSV Files (*.csv);;NumPy Files (*.npy);;All Files (*.*)"
        )
        
        if not file_path:
            return
        
        try:
            samples = self.lls_fitter.mcmc_results['samples']
            
            # Create a header for the columns
            header = "C0,C1,logNHI"
            
            # Check selected format
            if selected_filter == "CSV Files (*.csv)" or file_path.lower().endswith('.csv'):
                # Save as CSV
                np.savetxt(file_path, samples, delimiter=',', header=header)
                self.status_bar.showMessage(f"MCMC posterior saved to {file_path}", 3000)
                
            elif selected_filter == "NumPy Files (*.npy)" or file_path.lower().endswith('.npy'):
                # Save as NumPy array
                np.save(file_path, samples)
                self.status_bar.showMessage(f"MCMC posterior saved to {file_path}", 3000)
                
            else:
                # Default to CSV if extension is unknown
                np.savetxt(file_path, samples, delimiter=',', header=header)
                self.status_bar.showMessage(f"MCMC posterior saved to {file_path}", 3000)
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error saving MCMC posterior: {str(e)}")
    
    def save_results(self):
        """Save fit results to a JSON file"""
        if self.lls_fitter is None or (
            self.lls_fitter.curve_fit_results is None and 
            self.lls_fitter.mcmc_results is None
        ):
            QMessageBox.warning(self, "No Results", "No fit results available to save")
            return
        
        # Open file dialog for JSON
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Results", "", "JSON Files (*.json);;All Files (*.*)"
        )
        
        if not file_path:
            return
        
        try:
            import json
            import datetime
            
            # Get results and other data
            results = self.lls_fitter.get_results_summary()
            
            # Create simplified data structure
            data = {
                "metadata": {
                    "created": datetime.datetime.now().isoformat(),
                    "spectrum_file": self.file_path_edit.text(),
                    "redshift": float(self.redshift_edit.text())
                },
                "continuum_regions": self.get_current_continuum_regions(),
                "results": {}
            }
            
            # Add curve fit results if available
            if 'curve_fit' in results:
                data["results"]["curve_fit"] = results['curve_fit']
            
            # Add MCMC results if available
            if 'mcmc' in results:
                data["results"]["mcmc"] = results['mcmc']
                
                # Add percentiles if available
                if self.lls_fitter.mcmc_results and 'samples' in self.lls_fitter.mcmc_results:
                    samples = self.lls_fitter.mcmc_results['samples']
                    percentiles = {}
                    
                    for i, param_name in enumerate(['C0', 'C1', 'logNHI']):
                        p16, p50, p84 = np.percentile(samples[:, i], [16, 50, 84])
                        percentiles[param_name] = {
                            "p16": float(p16),
                            "p50": float(p50),
                            "p84": float(p84),
                            "upper_error": float(p84-p50),
                            "lower_error": float(p50-p16)
                        }
                    
                    data["results"]["mcmc_percentiles"] = percentiles
            
            # Write to file with nice formatting
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2)
            
            self.statusBar().showMessage(f"Results saved to {file_path}", 3000)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error saving results: {str(e)}")    

    def load_results(self):
        """Load and display saved results from a JSON file"""
        # Open file dialog
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Results", "", "JSON Files (*.json);;All Files (*.*)"
        )
        
        if not file_path:
            return
        
        try:
            import json
            
            # Read the JSON file
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Display a summary of the loaded results
            if 'metadata' in data and 'results' in data:
                # Create a dialog to show the loaded results
                dialog = QDialog(self)
                dialog.setWindowTitle("Loaded Results")
                dialog.resize(600, 400)
                
                layout = QVBoxLayout(dialog)
                
                # Create text display
                text_edit = QTextEdit()
                text_edit.setReadOnly(True)
                font = QFont("Monospace")
                font.setStyleHint(QFont.TypeWriter)
                text_edit.setFont(font)
                
                # Format results
                text = "Loaded LLS Fit Results\n"
                text += "====================\n\n"
                
                # Metadata
                meta = data['metadata']
                text += f"Spectrum: {meta.get('spectrum_file', 'N/A')}\n"
                text += f"Redshift: {meta.get('redshift', 'N/A')}\n"
                text += f"Created: {meta.get('created', 'N/A')}\n\n"
                
                # Continuum regions
                if 'continuum_regions' in data:
                    text += "Continuum Regions:\n"
                    for i, (wmin, wmax) in enumerate(data['continuum_regions']):
                        text += f"  Region {i+1}: {wmin:.1f} - {wmax:.1f} Å\n"
                    text += "\n"
                
                # Results
                if 'results' in data:
                    results = data['results']
                    
                    if 'curve_fit' in results:
                        cf = results['curve_fit']
                        text += "Curve Fit Results:\n"
                        text += f"  C0:      {cf['C0']:.6f} ± {cf['C0_err']:.6f}\n"
                        text += f"  C1:      {cf['C1']:.6f} ± {cf['C1_err']:.6f}\n"
                        text += f"  log N(HI): {cf['logNHI']:.4f} ± {cf['logNHI_err']:.4f}\n\n"
                    
                    if 'mcmc' in results:
                        mc = results['mcmc']
                        text += "MCMC Results:\n"
                        text += f"  C0:      {mc['C0']:.6f} ± {mc['C0_err']:.6f}\n"
                        text += f"  C1:      {mc['C1']:.6f} ± {mc['C1_err']:.6f}\n"
                        text += f"  log N(HI): {mc['logNHI']:.4f} ± {mc['logNHI_err']:.4f}\n\n"
                    
                    if 'mcmc_percentiles' in results:
                        percentiles = results['mcmc_percentiles']
                        text += "MCMC Percentiles:\n"
                        for param, values in percentiles.items():
                            p50, upper, lower = values['p50'], values['upper_error'], values['lower_error']
                            text += f"  {param}: {p50:.6f} (+{upper:.6f}) (-{lower:.6f})\n"
                
                text_edit.setText(text)
                layout.addWidget(text_edit)
                
                # Add buttons
                button_layout = QHBoxLayout()
                close_button = QPushButton("Close")
                close_button.clicked.connect(dialog.accept)
                copy_button = QPushButton("Copy to Clipboard")
                copy_button.clicked.connect(lambda: QApplication.clipboard().setText(text))
                
                button_layout.addWidget(copy_button)
                button_layout.addWidget(close_button)
                layout.addLayout(button_layout)
                
                # Show dialog
                dialog.exec_()
                
                self.statusBar().showMessage(f"Results loaded from {file_path}", 3000)
            else:
                QMessageBox.warning(self, "Invalid File", "The file does not contain valid LLS fit results")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading results: {str(e)}")
    
        
    def export_current_plot(self):
        """Export the current plot to an image file"""
        # Determine which plot is currently visible
        current_tab = None
        right_panel = None
        
        for widget in [self.preview_plot_widget, self.curve_fit_plot_widget, 
                      self.mcmc_plot_widget, self.corner_plot_widget]:
            parent = widget.parent()
            if isinstance(parent, QTabWidget) and parent.currentWidget() == widget:
                current_tab = widget
                right_panel = parent
                break
        
        if current_tab is None:
            QMessageBox.warning(self, "No Plot", "No plot is currently selected")
            return
        
        # Open file dialog
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Plot", "", "PNG Files (*.png);;PDF Files (*.pdf);;All Files (*.*)"
        )
        
        if not file_path:
            return
        
        try:
            # Save the figure
            current_tab.figure.savefig(file_path, dpi=300, bbox_inches='tight')
            self.statusBar().showMessage(f"Plot exported to {file_path}", 3000)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error exporting plot: {str(e)}")
    
    def view_detailed_results(self):
        """Show detailed results in a dialog"""
        if self.lls_fitter is None or (
            self.lls_fitter.curve_fit_results is None and 
            self.lls_fitter.mcmc_results is None
        ):
            QMessageBox.warning(self, "No Results", "No fit results available to view")
            return
        
        results = self.lls_fitter.get_results_summary()
        
        # Add samples to results if available
        if self.lls_fitter.mcmc_results and 'samples' in self.lls_fitter.mcmc_results:
            results['samples'] = self.lls_fitter.mcmc_results['samples']
        
        # Show dialog
        dialog = ResultsDialog(results, self)
        dialog.exec_()
    
    def show_about(self):
        """Show about dialog"""
        QMessageBox.about(self, "About LLS Fitter GUI",
                         """<b>LLS Fitter GUI</b><br><br>
                         A graphical interface for the LLSFitter class.<br><br>
                         This application provides tools for measuring the column density
                         of Lyman Limit Systems (LLS) in quasar spectra using both
                         curve-fitting and MCMC methods.""")
    
    def show_help(self):
        """Show help dialog"""
        help_text = """
        <h2>LLS Fitter GUI Help</h2>
        
        <h3>Basic Usage:</h3>
        <ol>
            <li>Load a spectrum using File → Open Spectrum or the Browse button</li>
            <li>Enter the absorption redshift</li>
            <li>Define continuum regions in the Continuum Regions tab</li>
            <li>Preview the continuum regions by clicking the Preview button</li>
            <li>Adjust fit parameters if needed in the Fit Parameters tab</li>
            <li>Run Curve Fit for a quick analysis or MCMC Fit for a more robust analysis</li>
            <li>View results in the plot tabs and in the results display at the bottom</li>
            <li>Export results or plots using the File menu</li>
        </ol>
        
        <h3>Continuum Regions:</h3>
        <p>Define wavelength regions where the continuum should be fitted. Regions should avoid
        absorption lines. The preview shows which points will be used in the fit.</p>
        
        <h3>Fit Parameters:</h3>
        <ul>
            <li><b>C0, C1:</b> Continuum parameters (constant and slope)</li>
            <li><b>log N(HI):</b> Logarithm of the neutral hydrogen column density</li>
            <li><b>Bounds:</b> Set minimum and maximum allowed values for each parameter</li>
            <li><b>MCMC Parameters:</b> Control the MCMC fitting process</li>
            <li><b>Plot Options:</b> Control how results are displayed</li>
        </ul>
        
        <h3>Results:</h3>
        <p>To view detailed results including parameter correlations, use Analysis → View Detailed Results
        or examine the corner plot tab after running MCMC.</p>
        """
        
        msg = QMessageBox(self)
        msg.setWindowTitle("LLS Fitter GUI Help")
        msg.setTextFormat(Qt.RichText)
        msg.setText(help_text)
        msg.exec_()
    
    

    def init_ui(self):
        # Create central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Create main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Create top controls
        top_controls = QWidget()
        top_layout = QHBoxLayout(top_controls)
        
        # File selection
        file_group = QGroupBox("Spectrum File")
        file_layout = QHBoxLayout(file_group)
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setReadOnly(True)
        browse_button = QPushButton("Browse...")
        browse_button.clicked.connect(self.browse_file)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(browse_button)
        
        ## Redshift input
        redshift_group = QGroupBox("Absorption Redshift")
        redshift_layout = QHBoxLayout(redshift_group)
        self.redshift_edit = QLineEdit()
        #self.redshift_edit.textChanged.connect(self.update_redshift)  # Add this line
        
        self.redshift_edit.returnPressed.connect(self.submit_redshift)
        redshift_submit_button = QPushButton("Submit")
        redshift_submit_button.clicked.connect(self.submit_redshift)
        redshift_layout.addWidget(self.redshift_edit)
        redshift_layout.addWidget(redshift_submit_button)
        
        # Add to top layout
        top_layout.addWidget(file_group, 3)
        top_layout.addWidget(redshift_group, 1)
        
        # Create splitter for main area
        splitter = QSplitter(Qt.Horizontal)
        
        # Left panel - controls
        left_panel = QTabWidget()
        
        # Tab 1: Continuum Regions
        continuum_tab = QWidget()
        continuum_layout = QVBoxLayout(continuum_tab)
        
        # Continuum regions table
        self.continuum_table = ContinuumRegionsTable()
        
        # Table controls
        table_controls = QHBoxLayout()
        add_button = QPushButton("Add Region")
        add_button.clicked.connect(self.continuum_table.add_row)
        remove_button = QPushButton("Remove Selected")
        remove_button.clicked.connect(self.continuum_table.remove_selected_rows)
        default_button = QPushButton("Default Regions")
        default_button.clicked.connect(self.continuum_table.set_default_regions)
        
        table_controls.addWidget(add_button)
        table_controls.addWidget(remove_button)
        table_controls.addWidget(default_button)
        
        continuum_layout.addWidget(self.continuum_table)
        continuum_layout.addLayout(table_controls)
        
        # Preview button
        preview_button = QPushButton("Preview Continuum Regions")
        preview_button.clicked.connect(self.preview_continuum_regions)
        continuum_layout.addWidget(preview_button)
        
        # Tab 2: Fit Parameters
        fit_tab = QWidget()
        fit_scroll = QScrollArea()
        fit_scroll.setWidgetResizable(True)
        fit_scroll_content = QWidget()
        fit_layout = QFormLayout(fit_scroll_content)
        fit_scroll.setWidget(fit_scroll_content)

        # Initial parameters
        params_group = QGroupBox("Initial Parameters")
        params_layout = QFormLayout(params_group)
        
        self.c0_input = QDoubleSpinBox()
        self.c0_input.setRange(-100, 100)
        self.c0_input.setValue(1.0)
        self.c0_input.setDecimals(4)
        
        self.c1_input = QDoubleSpinBox()
        self.c1_input.setRange(-10, 10)
        self.c1_input.setValue(0.0)
        self.c1_input.setDecimals(6)
        
        self.nhi_input = QDoubleSpinBox()
        self.nhi_input.setRange(14, 19)
        self.nhi_input.setValue(17.0)
        self.nhi_input.setDecimals(2)
        
        params_layout.addRow("C0:", self.c0_input)
        params_layout.addRow("C1:", self.c1_input)
        params_layout.addRow("log N(HI):", self.nhi_input)
        
        # Boundary conditions
        bounds_group = QGroupBox("Parameter Bounds")
        bounds_layout = QGridLayout(bounds_group)
        
        self.c0_min = QDoubleSpinBox()
        self.c0_min.setRange(-100, 100)
        self.c0_min.setValue(-10.0)
        self.c0_min.setDecimals(4)
        
        self.c0_max = QDoubleSpinBox()
        self.c0_max.setRange(-100, 100)
        self.c0_max.setValue(10.0)
        self.c0_max.setDecimals(4)
        
        self.c1_min = QDoubleSpinBox()
        self.c1_min.setRange(-10, 10)
        self.c1_min.setValue(-10.0)
        self.c1_min.setDecimals(6)
        
        self.c1_max = QDoubleSpinBox()
        self.c1_max.setRange(-10, 10)
        self.c1_max.setValue(10.0)
        self.c1_max.setDecimals(6)
        
        self.nhi_min = QDoubleSpinBox()
        self.nhi_min.setRange(13, 20)
        self.nhi_min.setValue(14.0)
        self.nhi_min.setDecimals(2)
        
        self.nhi_max = QDoubleSpinBox()
        self.nhi_max.setRange(13, 20)
        self.nhi_max.setValue(19.0)
        self.nhi_max.setDecimals(2)
        
        bounds_layout.addWidget(QLabel("Min"), 0, 1)
        bounds_layout.addWidget(QLabel("Max"), 0, 2)
        bounds_layout.addWidget(QLabel("C0:"), 1, 0)
        bounds_layout.addWidget(self.c0_min, 1, 1)
        bounds_layout.addWidget(self.c0_max, 1, 2)
        bounds_layout.addWidget(QLabel("C1:"), 2, 0)
        bounds_layout.addWidget(self.c1_min, 2, 1)
        bounds_layout.addWidget(self.c1_max, 2, 2)
        bounds_layout.addWidget(QLabel("log N(HI):"), 3, 0)
        bounds_layout.addWidget(self.nhi_min, 3, 1)
        bounds_layout.addWidget(self.nhi_max, 3, 2)
        
        # MCMC Parameters
        mcmc_group = QGroupBox("MCMC Parameters")
        mcmc_layout = QFormLayout(mcmc_group)
        
        self.walkers_input = QSpinBox()
        self.walkers_input.setRange(10, 500)
        self.walkers_input.setValue(50)
        
        self.steps_input = QSpinBox()
        self.steps_input.setRange(100, 10000)
        self.steps_input.setValue(500)
        
        self.burnin_input = QDoubleSpinBox()
        self.burnin_input.setRange(0.1, 0.5)
        self.burnin_input.setValue(0.2)
        self.burnin_input.setDecimals(2)
        self.burnin_input.setSingleStep(0.05)
        
        mcmc_layout.addRow("Number of Walkers:", self.walkers_input)
        mcmc_layout.addRow("Number of Steps:", self.steps_input)
        mcmc_layout.addRow("Burn-in Fraction:", self.burnin_input)
        
        # Plot Options
        plot_group = QGroupBox("Plot Options")
        plot_layout = QFormLayout(plot_group)  # Don't set the layout on the group box yet

        
        self.wmin_input = QDoubleSpinBox()
        self.wmin_input.setRange(800, 1000)
        self.wmin_input.setValue(880)
        
        self.wmax_input = QDoubleSpinBox()
        self.wmax_input.setRange(800, 1000)
        self.wmax_input.setValue(975)
        
        # Add y-axis limit inputs
        self.y_auto_checkbox = QCheckBox("Auto Y-Limits")
        self.y_auto_checkbox.setChecked(True)
        self.y_auto_checkbox.stateChanged.connect(self.toggle_y_limits)
        
        self.ymin_input = QDoubleSpinBox()
        self.ymin_input.setRange(-100, 100)
        self.ymin_input.setValue(0)
        self.ymin_input.setDecimals(1)
        self.ymin_input.setSingleStep(1)
        self.ymin_input.setEnabled(False)
        
        self.ymax_input = QDoubleSpinBox()
        self.ymax_input.setRange(0, 100)
        self.ymax_input.setValue(2)
        self.ymax_input.setDecimals(1)
        self.ymax_input.setSingleStep(1)
        self.ymax_input.setEnabled(False)
        
        self.show_regions_check = QCheckBox()
        self.show_regions_check.setChecked(True)
        
        self.show_realizations_check = QCheckBox()
        self.show_realizations_check.setChecked(True)
        
        self.n_realizations_input = QSpinBox()
        self.n_realizations_input.setRange(10, 500)
        self.n_realizations_input.setValue(100)
        
        # Add sigma spinbox
        self.sigma_input = QDoubleSpinBox()
        self.sigma_input.setRange(0.5, 10.0)
        self.sigma_input.setValue(3.0)
        self.sigma_input.setDecimals(1)
        self.sigma_input.setSingleStep(0.1)
        
        #Form rows
        plot_layout.addRow("Sigma Clip Threshold:", self.sigma_input)
        plot_layout.addRow("Min Wavelength (Å):", self.wmin_input)
        plot_layout.addRow("Max Wavelength (Å):", self.wmax_input)
        plot_layout.addRow("Auto Y-Limits:", self.y_auto_checkbox)
        plot_layout.addRow("Min Y Value:", self.ymin_input)
        plot_layout.addRow("Max Y Value:", self.ymax_input)
        plot_layout.addRow("Show Continuum Regions:", self.show_regions_check)
        plot_layout.addRow("Show MCMC Realizations:", self.show_realizations_check)
        plot_layout.addRow("Number of Realizations:", self.n_realizations_input)
             
        
        # Add groups to fit_layout
        fit_layout.addRow(params_group)
        fit_layout.addRow(bounds_group)
        fit_layout.addRow(mcmc_group)
        fit_layout.addRow(plot_group)

        fit_tab_layout = QVBoxLayout(fit_tab)
        fit_tab_layout.addWidget(fit_scroll)

        
        # Add tabs to the left panel
        left_panel.addTab(continuum_tab, "Continuum Regions")
        left_panel.addTab(fit_tab, "Fit Parameters")
        
        # Right panel - plots
        right_panel = QTabWidget()
        
        # Create plot widgets
        self.preview_plot_widget = MatplotlibWidget()
        self.curve_fit_plot_widget = MatplotlibWidget()
        self.mcmc_plot_widget = MatplotlibWidget()
        self.corner_plot_widget = MatplotlibWidget()
        
        # Add plot widgets to tabs
        right_panel.addTab(self.preview_plot_widget, "Preview")
        right_panel.addTab(self.curve_fit_plot_widget, "Curve Fit")
        right_panel.addTab(self.mcmc_plot_widget, "MCMC Fit")
        right_panel.addTab(self.corner_plot_widget, "Corner Plot")
        
        # Add panels to splitter
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([400, 800])  # Initial sizes
        
        # Bottom controls
        bottom_controls = QWidget()
        bottom_layout = QVBoxLayout(bottom_controls)
        
        # Button row
        button_row = QHBoxLayout()
        
        # Buttons for fitting
        self.curve_fit_button = QPushButton("Run Curve Fit")
        self.curve_fit_button.clicked.connect(self.run_curve_fit)
        self.curve_fit_button.setEnabled(False)
        
        self.mcmc_fit_button = QPushButton("Run MCMC Fit")
        self.mcmc_fit_button.clicked.connect(self.run_mcmc_fit)
        self.mcmc_fit_button.setEnabled(False)
        
        # Results display
        self.result_label = QLabel("No fit results yet")
        self.result_label.setAlignment(Qt.AlignCenter)
        self.result_label.setFrameShape(QFrame.Panel)
        self.result_label.setFrameShadow(QFrame.Sunken)
        font = QFont()
        font.setBold(True)
        self.result_label.setFont(font)
        
        # Add to button row
        button_row.addWidget(self.curve_fit_button)
        button_row.addWidget(self.mcmc_fit_button)
        button_row.addWidget(self.result_label)
        
        # Progress bar
        #self.progress_bar = QProgressBar()
        #self.progress_bar.setRange(0, 100)
        #self.progress_bar.setValue(0)
        #self.progress_bar.setVisible(False)
        
        # Add to bottom layout
        bottom_layout.addLayout(button_row)
        #bottom_layout.addWidget(self.progress_bar)
        
        # Add all widgets to main layout
        main_layout.addWidget(top_controls)
        main_layout.addWidget(splitter)
        main_layout.addWidget(bottom_controls)
        
        # Set up status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        
        # Set initial state
        self.continuum_table.set_default_regions()

    def submit_redshift(self):
        """Update the redshift when the submit button is clicked or Enter is pressed"""
        if self.lls_fitter is not None:
            try:
                zabs_text = self.redshift_edit.text()
                if zabs_text:
                    zabs = float(zabs_text)
                    self.lls_fitter.set_redshift(zabs)
                    self.status_bar.showMessage(f"Redshift updated to {zabs:.6f}", 3000)
                    self.preview_continuum_regions()
                else:
                    self.status_bar.showMessage("Please enter a valid redshift", 3000)
            except ValueError:
                self.status_bar.showMessage("Invalid redshift format", 3000)
                QMessageBox.warning(self, "Invalid Input", "Please enter a valid numeric redshift value")

    def toggle_y_limits(self, state):
        """Enable or disable custom y-limit inputs based on checkbox state"""
        enabled = not bool(state)  # Checkbox checked = auto limits = disable inputs
        self.ymin_input.setEnabled(enabled)
        self.ymax_input.setEnabled(enabled)
            
    def browse_file(self):
        """Open a file dialog to select a spectrum file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Spectrum File", "", "FITS Files (*.fits);;All Files (*.*)"
        )
        
        if file_path:
            self.file_path_edit.setText(file_path)
            self.load_spectrum()
    
    def load_spectrum(self):
        """Load the spectrum from the selected file"""
        file_path = self.file_path_edit.text()
        if not file_path or not os.path.exists(file_path):
            QMessageBox.warning(self, "Error", "Please select a valid spectrum file")
            return
        
        try:
            zabs_text = self.redshift_edit.text()
            if not zabs_text:
                QMessageBox.warning(self, "Error", "Please enter an absorption redshift")
                return
            
            zabs = float(zabs_text)
            
            # Create LLSFitter object
            self.lls_fitter = LLSFitter(file_path, zabs)
            
            # Update status
            self.statusBar().showMessage(f"Loaded spectrum from {file_path} with z_abs = {zabs:.6f}")
            
            # Enable controls
            self.curve_fit_button.setEnabled(True)
            self.mcmc_fit_button.setEnabled(True)
            
            # Preview continuum regions
            self.preview_continuum_regions()
            
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid numeric redshift value")
            return
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading spectrum: {str(e)}")
            self.statusBar().showMessage("Error loading spectrum")

    
    def get_current_continuum_regions(self):
        """Get the current continuum regions from the table"""
        return self.continuum_table.get_regions()
    
    def update_lls_fitter_params(self):
        """Update the LLSFitter object with the current GUI parameters"""
        if self.lls_fitter is None:
            return False
        
        try:
            # Set continuum regions
            regions = self.get_current_continuum_regions()
            if not regions:
                QMessageBox.warning(self, "Error", "No valid continuum regions defined")
                return False
            
            self.lls_fitter.set_continuum_regions(regions)
            
            # Set initial parameters
            theta_init = np.array([
                self.c0_input.value(),
                self.c1_input.value(),
                self.nhi_input.value()
            ])
            
            # Set bounds
            self.lls_fitter.bounds_lower = np.array([
                self.c0_min.value(),
                self.c1_min.value(),
                self.nhi_min.value()
            ])
            
            self.lls_fitter.bounds_upper = np.array([
                self.c0_max.value(),
                self.c1_max.value(),
                self.nhi_max.value()
            ])
            
            # Set the initial parameters
            self.lls_fitter.theta_init = theta_init

            # Set sigma clip value
            self.lls_fitter.sigma_clip = self.sigma_input.value()
        
            
            return True
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error updating parameters: {str(e)}")
            return False
    
    def preview_continuum_regions(self):
        """Preview the continuum regions on the spectrum"""
        if self.lls_fitter is None:
            self.load_spectrum()
            if self.lls_fitter is None:  # If still None, loading failed
                return
        
        try:
            # Update the LLSFitter with current parameters
            if not self.update_lls_fitter_params():
                return
            
            # Capture printed warnings by temporarily redirecting stdout
            import io
            import sys
            original_stdout = sys.stdout
            captured_output = io.StringIO()
            sys.stdout = captured_output
            
            # Get the continuum mask
            mask = self.lls_fitter.get_continuum_mask()
            
            # Restore stdout
            sys.stdout = original_stdout
            warnings = captured_output.getvalue()
            
            # Check for warnings and show in message box if present
            if warnings:
                QMessageBox.warning(self, "Continuum Region Warnings", warnings)
            
            # Check if mask is empty (no valid continuum points)
            if not np.any(mask):
                QMessageBox.critical(self, "Error", "No valid continuum points found. Cannot proceed with fitting.")
                return
            
            # Clear the figure
            self.preview_plot_widget.clear_figure()
            
            # Get figure and add new axis
            fig = self.preview_plot_widget.get_figure()
            ax = fig.add_subplot(111)
            
            # Get the wavelength range
            wmin = self.wmin_input.value()
            wmax = self.wmax_input.value()
            
            # Get y-limits if not using auto
            ymin = None
            ymax = None
            if not self.y_auto_checkbox.isChecked():
                ymin = self.ymin_input.value()
                ymax = self.ymax_input.value()
            
            # Get the continuum mask
            mask = self.lls_fitter.get_continuum_mask()
            
            # Plot spectrum
            ax.step(self.lls_fitter.rest_wave, self.lls_fitter.flux, 'b-', 
                    alpha=0.7, where='mid', label='Spectrum')
            
            # Plot continuum points
            ax.plot(self.lls_fitter.rest_wave[mask], self.lls_fitter.flux[mask], 'r.', 
                   alpha=0.7, markersize=3, label='Continuum Points')
            
            # Add vertical line at Lyman limit
            ax.axvline(x=912, color='k', linestyle='--', alpha=0.7, label='Lyman Limit (912 Å)')
            
            # Highlight continuum regions
            if self.show_regions_check.isChecked():
                y_box_min = ymin if ymin is not None else 0
                y_box_max = ymax if ymax is not None else 2  # Default, will be adjusted later
                
                for rmin, rmax in self.get_current_continuum_regions():
                    rect = plt.Rectangle((rmin, y_box_min), rmax-rmin, y_box_max-y_box_min, 
                                        color='gray', alpha=0.15, zorder=0)
                    ax.add_patch(rect)
            
            # Set x-axis limits
            ax.set_xlim(wmin, wmax)
            
            # Set y-axis limits (either custom or auto-calculated)
            if ymin is not None and ymax is not None:
                ax.set_ylim(ymin, ymax)
                
                # Adjust rectangle heights for continuum regions
                if self.show_regions_check.isChecked():
                    for patch in ax.patches:
                        if isinstance(patch, plt.Rectangle):
                            patch.set_height(ymax - ymin)
                            patch.set_y(ymin)
            else:
                # Auto-calculate y limits
                visible_flux = self.lls_fitter.flux[
                    (self.lls_fitter.rest_wave >= wmin) & (self.lls_fitter.rest_wave <= wmax)
                ]
                auto_ymin = np.percentile(visible_flux, 1)
                auto_ymax = np.percentile(visible_flux, 99)
                margin = 0.2 * (auto_ymax - auto_ymin)
                ax.set_ylim(max(0, auto_ymin - margin), auto_ymax + margin)
                
                # Adjust rectangle heights for continuum regions
                if self.show_regions_check.isChecked():
                    for patch in ax.patches:
                        if isinstance(patch, plt.Rectangle):
                            patch.set_height(ax.get_ylim()[1] - ax.get_ylim()[0])
                            patch.set_y(ax.get_ylim()[0])
            
            # Add labels and legend
            ax.set_xlabel('Rest Wavelength (Å)')
            ax.set_ylabel('Normalized Flux')
            ax.set_title('Continuum Regions Preview')
            ax.legend(loc='best')
            ax.grid(True, alpha=0.3)
            
            # Update the canvas
            self.preview_plot_widget.canvas.draw()
            
            # Switch to the preview tab
            right_panel = self.preview_plot_widget.parent()
            if isinstance(right_panel, QTabWidget):
                right_panel.setCurrentWidget(self.preview_plot_widget)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error previewing continuum regions: {str(e)}")
    
    def run_curve_fit(self):
        """Run the curve_fit process"""
        if self.lls_fitter is None:
            QMessageBox.warning(self, "Error", "No spectrum loaded")
            return
        
        try:
            # Update the LLSFitter with current parameters
            if not self.update_lls_fitter_params():
                return
            
            # Capture printed warnings by temporarily redirecting stdout
            import io
            import sys
            original_stdout = sys.stdout
            captured_output = io.StringIO()
            sys.stdout = captured_output
            
            # Get the continuum mask
            mask = self.lls_fitter.get_continuum_mask()
            
            # Restore stdout
            sys.stdout = original_stdout
            warnings = captured_output.getvalue()
            
            # Check for warnings and show in message box if present
            if warnings:
                QMessageBox.warning(self, "Continuum Region Warnings", warnings)
            
            # Check if mask is empty (no valid continuum points)
            if not np.any(mask):
                QMessageBox.critical(self, "Error", "No valid continuum points found. Cannot proceed with fitting.")
                return

            # Show progress in status bar
            self.status_bar.showMessage("Running curve_fit...")
            #self.progress_bar.setVisible(True)
            #self.progress_bar.setValue(20)  # Initial progress
            QApplication.processEvents()  # Update UI
            
            # Disable UI elements during fitting
            self.curve_fit_button.setEnabled(False)
            self.mcmc_fit_button.setEnabled(False)
            
            # Run curve_fit
            popt, pcov = self.lls_fitter.fit_curve_fit(self.lls_fitter.theta_init)
            
            # Update progress
            #self.progress_bar.setValue(70)
            QApplication.processEvents()

            # Get plot limits
            wmin = self.wmin_input.value()
            wmax = self.wmax_input.value()
            ymin = None
            ymax = None
            if not self.y_auto_checkbox.isChecked():
                ymin = self.ymin_input.value()
                ymax = self.ymax_input.value()
            
            # Plot the results
            self.curve_fit_plot_widget.clear_figure()
            fig, ax = self.lls_fitter.plot_fit(
                method='curve_fit',
                show_continuum_regions=self.show_regions_check.isChecked(),
                wmin=wmin,
                wmax=wmax,
                ymin=ymin,
                ymax=ymax,
                figsize=(8, 6)
            )
            
            # Replace the matplotlib figure in the widget
            self.curve_fit_plot_widget.figure = fig
            self.curve_fit_plot_widget.canvas.figure = fig
            self.curve_fit_plot_widget.canvas.draw()
            
            # Update progress
            #self.progress_bar.setValue(100)
            QApplication.processEvents()
            
            # Switch to the curve_fit tab
            right_panel = self.curve_fit_plot_widget.parent()
            if isinstance(right_panel, QTabWidget):
                right_panel.setCurrentWidget(self.curve_fit_plot_widget)
            
            # Update results display
            logNHI = popt[2]
            logNHI_err = np.sqrt(pcov[2, 2])
            self.result_label.setText(f"Curve Fit Results: log N(HI) = {logNHI:.2f} ± {logNHI_err:.2f}")
            
            # Update status
            self.status_bar.showMessage("Curve fit completed successfully", 3000)
            
            # Hide progress bar after a delay
            #QTimer.singleShot(1000, lambda: self.progress_bar.setVisible(False))
            
            # Update initial parameters for potential MCMC run
            self.c0_input.setValue(popt[0])
            self.c1_input.setValue(popt[1])
            self.nhi_input.setValue(popt[2])
            
            # Re-enable UI elements
            self.curve_fit_button.setEnabled(True)
            self.mcmc_fit_button.setEnabled(True)
            
        except Exception as e:
            # Handle errors
            self.status_bar.showMessage("Error in curve_fit")
            #self.progress_bar.setVisible(False)
            self.curve_fit_button.setEnabled(True)
            self.mcmc_fit_button.setEnabled(True)
            QMessageBox.critical(self, "Error", f"Error in curve_fit: {str(e)}")
    
    def run_mcmc_fit(self):
        """Run the MCMC fitting process"""
        if self.lls_fitter is None:
            QMessageBox.warning(self, "Error", "No spectrum loaded")
            return
        
        try:
            # Update the LLSFitter with current parameters
            if not self.update_lls_fitter_params():
                return

            # Capture printed warnings by temporarily redirecting stdout
            import io
            import sys
            original_stdout = sys.stdout
            captured_output = io.StringIO()
            sys.stdout = captured_output
            
            # Get the continuum mask
            mask = self.lls_fitter.get_continuum_mask()
            
            # Restore stdout
            sys.stdout = original_stdout
            warnings = captured_output.getvalue()
            
            # Check for warnings and show in message box if present
            if warnings:
                QMessageBox.warning(self, "Continuum Region Warnings", warnings)
            
            # Check if mask is empty (no valid continuum points)
            if not np.any(mask):
                QMessageBox.critical(self, "Error", "No valid continuum points found. Cannot proceed with fitting.")
                return

            
            # Confirm for long runs
            nwalkers = self.walkers_input.value()
            nsteps = self.steps_input.value()
            
            if nwalkers * nsteps > 50000:
                result = QMessageBox.question(
                    self, "Confirm MCMC Run",
                    f"You are about to run MCMC with {nwalkers} walkers and {nsteps} steps, "
                    f"which may take a long time. Continue?",
                    QMessageBox.Yes | QMessageBox.No
                )
                
                if result == QMessageBox.No:
                    return
            
            # Show progress in status bar
            self.status_bar.showMessage("Running MCMC fit (this may take a while)...")
            #self.progress_bar.setVisible(True)
            #self.progress_bar.setValue(0)  # Initial progress
            QApplication.processEvents()  # Update UI
            
            # Disable UI elements during fitting
            self.curve_fit_button.setEnabled(False)
            self.mcmc_fit_button.setEnabled(False)
            
            # Define a custom progress callback for MCMC
            # This will be called by the modified LLSFitter.fit_emcee method
            #def progress_callback(iteration, total_iterations):
            #    progress = int(100 * iteration / total_iterations)
            #    self.progress_bar.setValue(progress)
            #    if iteration % 10 == 0:  # Only update UI every 10 iterations
            #        self.status_bar.showMessage(f"MCMC progress: {progress}% ({iteration}/{total_iterations} steps)")
            #        QApplication.processEvents()
            
            # Run MCMC
            # Note: We'll need to modify LLSFitter.fit_emcee to accept a progress_callback
            sampler, samples = self.lls_fitter.fit_emcee(
                nwalkers=nwalkers,
                nsteps=nsteps,
                burnin_frac=self.burnin_input.value(),
                theta_init=self.lls_fitter.theta_init
                # progress_callback=progress_callback  # Uncomment after modifying LLSFitter
            )
            
            # Since we don't have the callback in the original class, simulate progress updates
            # Remove this block after implementing the callback in LLSFitter
            for i in range(0, 101, 10):
                #self.progress_bar.setValue(i)
                self.status_bar.showMessage(f"MCMC progress: {i}%")
                QApplication.processEvents()
                QTimer.singleShot(100, lambda: None)  # Small delay for visual feedback
            
            # Update progress
            #self.progress_bar.setValue(80)
            self.status_bar.showMessage("Generating plots...")
            QApplication.processEvents()


            # Get plot limits
            wmin = self.wmin_input.value()
            wmax = self.wmax_input.value()
            ymin = None
            ymax = None
            if not self.y_auto_checkbox.isChecked():
                ymin = self.ymin_input.value()
                ymax = self.ymax_input.value()

            
            # Plot the MCMC fit
            self.mcmc_plot_widget.clear_figure()
            fig, ax = self.lls_fitter.plot_fit(
                method='mcmc',
                show_continuum_regions=self.show_regions_check.isChecked(),
                wmin=wmin,
                wmax=wmax,
                ymin=ymin,
                ymax=ymax,
                show_realizations=self.show_realizations_check.isChecked(),
                n_realizations=self.n_realizations_input.value(),
                figsize=(8, 6)
            )
            
            # Replace the matplotlib figure in the widget
            self.mcmc_plot_widget.figure = fig
            self.mcmc_plot_widget.canvas.figure = fig
            self.mcmc_plot_widget.canvas.draw()
            
            # Update progress
            #self.progress_bar.setValue(90)
            self.status_bar.showMessage("Generating corner plot...")
            QApplication.processEvents()
            
            # Plot the corner plot
            self.corner_plot_widget.clear_figure()
            corner_fig = self.lls_fitter.plot_corner(figsize=(8, 8))
            
            # Replace the matplotlib figure in the widget
            self.corner_plot_widget.figure = corner_fig
            self.corner_plot_widget.canvas.figure = corner_fig
            self.corner_plot_widget.canvas.draw()
            
            # Update progress
            #self.progress_bar.setValue(100)
            QApplication.processEvents()
            
            # Switch to the MCMC tab
            right_panel = self.mcmc_plot_widget.parent()
            if isinstance(right_panel, QTabWidget):
                right_panel.setCurrentWidget(self.mcmc_plot_widget)
            
            # Get results
            results = self.lls_fitter.get_results_summary()
            logNHI = results['mcmc']['logNHI']
            logNHI_err = results['mcmc']['logNHI_err']
            
            # Update results display
            self.result_label.setText(f"MCMC Results: log N(HI) = {logNHI:.2f} ± {logNHI_err:.2f}")
            
            # Update status
            self.status_bar.showMessage("MCMC fit completed successfully", 3000)
            
            # Hide progress bar after a delay
            #QTimer.singleShot(1000, lambda: self.progress_bar.setVisible(False))
            
            # Re-enable UI elements
            self.curve_fit_button.setEnabled(True)
            self.mcmc_fit_button.setEnabled(True)
            
        except Exception as e:
            # Handle errors
            self.status_bar.showMessage("Error in MCMC fit")
            #self.progress_bar.setVisible(False)
            self.curve_fit_button.setEnabled(True)
            self.mcmc_fit_button.setEnabled(True)
            QMessageBox.critical(self, "Error", f"Error in MCMC fit: {str(e)}")
    
    def save_settings(self):
        """Save current settings"""
        settings = self.settings
        
        # Save redshift
        settings.setValue("redshift", self.redshift_edit.text())
        
        # Save fit parameters
        settings.setValue("c0", self.c0_input.value())
        settings.setValue("c1", self.c1_input.value())
        settings.setValue("nhi", self.nhi_input.value())
        
        # Save bounds
        settings.setValue("c0_min", self.c0_min.value())
        settings.setValue("c0_max", self.c0_max.value())
        settings.setValue("c1_min", self.c1_min.value())
        settings.setValue("c1_max", self.c1_max.value())
        settings.setValue("nhi_min", self.nhi_min.value())
        settings.setValue("nhi_max", self.nhi_max.value())
        
        # Save MCMC parameters
        settings.setValue("walkers", self.walkers_input.value())
        settings.setValue("steps", self.steps_input.value())
        settings.setValue("burnin", self.burnin_input.value())
        
        # Save plot options
        settings.setValue("wmin", self.wmin_input.value())
        settings.setValue("wmax", self.wmax_input.value())
        settings.setValue("y_auto", self.y_auto_checkbox.isChecked())
        settings.setValue("ymin", self.ymin_input.value())
        settings.setValue("ymax", self.ymax_input.value())
        settings.setValue("show_regions", self.show_regions_check.isChecked())
        settings.setValue("show_realizations", self.show_realizations_check.isChecked())
        settings.setValue("n_realizations", self.n_realizations_input.value())
    
        
        # Save window geometry
        settings.setValue("geometry", self.saveGeometry())
        settings.setValue("state", self.saveState())

        # Save sigma value
        settings.setValue("sigma_clip", self.sigma_input.value())
    
    def load_settings(self):
        """Load saved settings"""
        settings = self.settings
        
        # Load redshift
        self.redshift_edit.setText(settings.value("redshift", ""))
        
        # Load fit parameters
        self.c0_input.setValue(float(settings.value("c0", 1.0)))
        self.c1_input.setValue(float(settings.value("c1", 0.0)))
        self.nhi_input.setValue(float(settings.value("nhi", 17.0)))
        
        # Load bounds
        self.c0_min.setValue(float(settings.value("c0_min", -10.0)))
        self.c0_max.setValue(float(settings.value("c0_max", 10.0)))
        self.c1_min.setValue(float(settings.value("c1_min", -10.0)))
        self.c1_max.setValue(float(settings.value("c1_max", 10.0)))
        self.nhi_min.setValue(float(settings.value("nhi_min", 14.0)))
        self.nhi_max.setValue(float(settings.value("nhi_max", 20.0)))
        
        # Load MCMC parameters
        self.walkers_input.setValue(int(settings.value("walkers", 50)))
        self.steps_input.setValue(int(settings.value("steps", 500)))
        self.burnin_input.setValue(float(settings.value("burnin", 0.2)))
        
        # Load plot options
        self.wmin_input.setValue(float(settings.value("wmin", 880)))
        self.wmax_input.setValue(float(settings.value("wmax", 975)))
        y_auto = settings.value("y_auto", "true") == "true"
        self.y_auto_checkbox.setChecked(y_auto)
        self.ymin_input.setValue(float(settings.value("ymin", 0.0)))
        self.ymax_input.setValue(float(settings.value("ymax", 2.0)))
        self.ymin_input.setEnabled(not y_auto)
        self.ymax_input.setEnabled(not y_auto)
        self.show_regions_check.setChecked(settings.value("show_regions", "true") == "true")
        self.show_realizations_check.setChecked(settings.value("show_realizations", "true") == "true")
        self.n_realizations_input.setValue(int(settings.value("n_realizations", 100)))
        
        # Load sigma value
        self.sigma_input.setValue(float(settings.value("sigma_clip", 3.0)))


        # Load window geometry
        geometry = settings.value("geometry")
        if geometry:
            self.restoreGeometry(geometry)
        state = settings.value("state")
        if state:
            self.restoreState(state)
    
    def closeEvent(self, event):
        """Handle close event to save settings"""
        self.save_settings()
        super(LLSFitterGUI, self).closeEvent(event)


def main():
    app = QApplication(sys.argv)
    window = LLSFitterGUI()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
