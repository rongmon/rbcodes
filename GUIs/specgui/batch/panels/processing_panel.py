# rbcodes/GUIs/specgui/batch/panels/processing_panel.py
import os
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                          QPushButton, QProgressBar, QGroupBox, QTableWidget,
                          QTableWidgetItem, QHeaderView, QCheckBox, QSpinBox,
                          QComboBox, QFormLayout, QScrollArea, QTextEdit,
                          QMessageBox, QFileDialog)
from PyQt5.QtCore import pyqtSignal, Qt, QThread, QObject
from PyQt5.QtGui import QColor


class BatchProcessor(QObject):
    """Worker object that processes batch items in a separate thread."""
    
    # Signals
    progress = pyqtSignal(int, int)  # current, total
    item_completed = pyqtSignal(int, dict)  # index, results
    item_failed = pyqtSignal(int, str)  # index, error
    finished = pyqtSignal(bool, list)  # success, results
    
    def __init__(self, controller, selected_indices=None):
        super().__init__()
        self.controller = controller
        self.selected_indices = selected_indices
        self.cancelled = False
    
    def process(self):
        """Process batch items."""
        success, results = self.controller.process_batch(self.selected_indices)
        self.finished.emit(success, results)
    
    def cancel(self):
        """Cancel batch processing."""
        self.cancelled = True


class ProcessingPanel(QWidget):
    """Panel for configuring and running batch processing."""
    
    # Signals
    batch_started = pyqtSignal()  # Emitted when batch processing starts
    batch_completed = pyqtSignal(bool, list)  # success, results
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.selected_indices = None  # Indices of items to process
        self.processor_thread = None  # Thread for batch processing
        self.processor = None  # Batch processor object
        self.init_ui()
        
        # Connect to controller signals
        self.controller.batch_progress.connect(self.update_progress)
        self.controller.batch_item_completed.connect(self.on_item_completed)
        self.controller.batch_item_failed.connect(self.on_item_failed)
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Processing settings
        settings_group = QGroupBox("Processing Settings")
        settings_layout = QFormLayout()
        
        # Continuum fitting method
        self.cont_method = QComboBox()
        self.cont_method.addItems(["Polynomial", "Flat (Skip Fitting)", "Interactive"])
        self.cont_method.currentTextChanged.connect(self.on_cont_method_changed)
        settings_layout.addRow("Continuum Method:", self.cont_method)
        
        # Polynomial order
        self.poly_order = QSpinBox()
        self.poly_order.setRange(1, 10)
        self.poly_order.setValue(3)
        settings_layout.addRow("Polynomial Order:", self.poly_order)
        
        # Use weights
        self.use_weights = QCheckBox()
        settings_layout.addRow("Use Weights:", self.use_weights)
        
        # Calculate SNR
        self.calc_snr = QCheckBox()
        self.calc_snr.setChecked(True)
        settings_layout.addRow("Calculate SNR:", self.calc_snr)
        
        # Binsize for SNR
        self.binsize = QSpinBox()
        self.binsize.setRange(1, 10)
        self.binsize.setValue(3)
        settings_layout.addRow("SNR Bin Size:", self.binsize)
        
        # Save individual JSON files
        self.save_json = QCheckBox()
        self.save_json.setChecked(True)
        settings_layout.addRow("Save Individual JSON Files:", self.save_json)
        
        # Output directory
        dir_layout = QHBoxLayout()
        self.output_dir = QLabel("Not set")
        self.browse_dir_btn = QPushButton("Browse")
        self.browse_dir_btn.clicked.connect(self.browse_output_dir)
        
        dir_layout.addWidget(self.output_dir)
        dir_layout.addWidget(self.browse_dir_btn)
        
        settings_layout.addRow("Output Directory:", dir_layout)
        
        # Apply settings button
        self.apply_settings_btn = QPushButton("Apply Settings")
        self.apply_settings_btn.clicked.connect(self.apply_settings)
        settings_layout.addRow("", self.apply_settings_btn)
        
        settings_group.setLayout(settings_layout)
        main_layout.addWidget(settings_group)
        
        # Progress section
        progress_group = QGroupBox("Batch Progress")
        progress_layout = QVBoxLayout()
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)
        
        # Status label
        self.status_label = QLabel("Ready to process batch items.")
        progress_layout.addWidget(self.status_label)
        
        # Processing buttons
        buttons_layout = QHBoxLayout()
        
        self.run_btn = QPushButton("Run Batch Process")
        self.run_btn.clicked.connect(self.run_batch)
        
        self.run_selected_btn = QPushButton("Run Selected Items")
        self.run_selected_btn.clicked.connect(self.run_selected)
        self.run_selected_btn.setEnabled(False)  # Initially disabled
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.cancel_processing)
        self.cancel_btn.setEnabled(False)  # Initially disabled
        
        buttons_layout.addWidget(self.run_btn)
        buttons_layout.addWidget(self.run_selected_btn)
        buttons_layout.addWidget(self.cancel_btn)
        
        progress_layout.addLayout(buttons_layout)
        
        progress_group.setLayout(progress_layout)
        main_layout.addWidget(progress_group)
        
        # Error log
        error_group = QGroupBox("Error Log")
        error_layout = QVBoxLayout()
        
        self.error_log = QTextEdit()
        self.error_log.setReadOnly(True)
        error_layout.addWidget(self.error_log)
        
        error_group.setLayout(error_layout)
        main_layout.addWidget(error_group)
    
    def on_cont_method_changed(self, method):
        """Handle changes to the continuum fitting method."""
        # Enable/disable polynomial order based on method
        is_polynomial = method == "Polynomial"
        self.poly_order.setEnabled(is_polynomial)
    
    def browse_output_dir(self):
        """Browse for output directory."""
        directory = QFileDialog.getExistingDirectory(
            self, "Select Output Directory", "", 
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks
        )
        
        if directory:
            self.output_dir.setText(directory)
    
    def apply_settings(self):
        """Apply the current settings to the controller."""
        # Get method
        method = self.cont_method.currentText().lower()
        if "polynomial" in method:
            method = "polynomial"
        elif "flat" in method:
            method = "flat"
        else:
            method = "interactive"
        
        # Update settings
        settings = {
            'continuum_method': method,
            'polynomial_order': self.poly_order.value(),
            'use_weights': self.use_weights.isChecked(),
            'calculate_snr': self.calc_snr.isChecked(),
            'binsize': self.binsize.value(),
            'save_individual_json': self.save_json.isChecked(),
            'output_directory': self.output_dir.text() if self.output_dir.text() != "Not set" else ""
        }
        
        self.controller.update_batch_settings(settings)
        
        QMessageBox.information(
            self, "Settings Applied", 
            "Batch processing settings have been applied."
        )
    
    def set_selected_items(self, indices):
        """Set the indices of selected items to process."""
        self.selected_indices = indices
        self.run_selected_btn.setEnabled(len(indices) > 0)
    
    def run_batch(self):
        """Run batch processing on all items."""
        self._start_processing()
    
    def run_selected(self):
        """Run batch processing on selected items."""
        if not self.selected_indices:
            QMessageBox.warning(
                self, "No Selection", 
                "Please select items to process in the Review tab."
            )
            return
        
        self._start_processing(self.selected_indices)
    
    def _start_processing(self, selected_indices=None):
        """Start the batch processing thread."""
        # Check if output directory is set if saving JSON files
        if (self.save_json.isChecked() and 
            (self.output_dir.text() == "Not set" or not os.path.isdir(self.output_dir.text()))):
            QMessageBox.warning(
                self, "Output Directory Not Set", 
                "Please set a valid output directory for saving JSON files."
            )
            return
        
        # Apply current settings
        self.apply_settings()
        
        # Create processor and thread
        self.processor = BatchProcessor(self.controller, selected_indices)
        self.processor_thread = QThread()
        self.processor.moveToThread(self.processor_thread)
        
        # Connect signals
        self.processor_thread.started.connect(self.processor.process)
        self.processor.progress.connect(self.update_progress)
        self.processor.item_completed.connect(self.on_item_completed)
        self.processor.item_failed.connect(self.on_item_failed)
        self.processor.finished.connect(self._on_processing_finished)
        
        # Update UI
        self.run_btn.setEnabled(False)
        self.run_selected_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.progress_bar.setValue(0)
        self.error_log.clear()
        
        # Emit signal that batch processing started
        self.batch_started.emit()
        
        # Start thread
        self.processor_thread.start()
    
    def _on_processing_finished(self, success, results):
        """Handle processing finished signal."""
        # Clean up thread
        self.processor_thread.quit()
        self.processor_thread.wait()
        
        # Update UI
        self.run_btn.setEnabled(True)
        self.run_selected_btn.setEnabled(self.selected_indices is not None and len(self.selected_indices) > 0)
        self.cancel_btn.setEnabled(False)
        
        # Emit signal that batch processing completed
        self.batch_completed.emit(success, results)
    
    def cancel_processing(self):
        """Cancel the current batch processing."""
        if self.processor:
            self.processor.cancel()
            self.status_label.setText("Cancelling batch processing...")
    
    def update_progress(self, current, total):
        """Update the progress bar."""
        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(current)
        self.status_label.setText(f"Processing item {current} of {total}...")
    
    def on_item_completed(self, index, result):
        """Handle item completed signal."""
        # Could update a results table here
        pass
    
    def on_item_failed(self, index, error):
        """Handle item failed signal."""
        # Add to error log
        self.error_log.append(f"Error processing item {index+1}: {error}")