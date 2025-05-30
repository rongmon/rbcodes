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
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.cancelled = False
    
    def process(self):
        """Process batch items."""
        success, results = self.controller.process_batch(None)  # Always process all items
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
        
        # Instructions
        instructions = QLabel(
            "Configure processing settings that will be applied to all batch items. "
            "Click 'Apply Settings' to save your choices, then 'Run Batch Process' to start."
        )
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)
        
        # Processing settings
        settings_group = QGroupBox("Processing Settings")
        settings_layout = QFormLayout()
        
        # Continuum fitting method
        self.cont_method = QComboBox()
        self.cont_method.addItems(["Polynomial", "Flat (Skip Fitting)"])
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
        
        # Apply Settings button
        self.apply_settings_btn = QPushButton("Apply Settings to All Items")
        self.apply_settings_btn.clicked.connect(self.apply_settings_to_master_table)
        self.apply_settings_btn.setToolTip("Save these settings and apply them to all batch items")
        settings_layout.addRow("", self.apply_settings_btn)
        
        # Status for settings
        self.settings_status = QLabel("Settings not applied yet.")
        self.settings_status.setStyleSheet("QLabel { color: red; font-style: italic; }")
        settings_layout.addRow("Status:", self.settings_status)
        
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
        
        # Processing button (compact size)
        self.run_btn = QPushButton("Run Batch Process")
        self.run_btn.clicked.connect(self.run_batch)
        self.run_btn.setMaximumWidth(200)  # Make button smaller
        
        # Center the button
        button_container = QHBoxLayout()
        button_container.addStretch()
        button_container.addWidget(self.run_btn)
        button_container.addStretch()
        
        progress_layout.addLayout(button_container)
        
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
        
        # Mark settings as not applied
        self.settings_status.setText("Settings changed - click 'Apply Settings' to save.")
        self.settings_status.setStyleSheet("QLabel { color: red; font-style: italic; }")
    
    def apply_settings_to_master_table(self):
        """Apply current UI settings to all items in the master table."""
        try:
            # Get settings from UI
            method = self.cont_method.currentText().lower()
            if "polynomial" in method:
                continuum_method = "polynomial"
            else:
                continuum_method = "flat"
            
            settings = {
                'continuum_method': continuum_method,
                'continuum_order': self.poly_order.value(),
                'use_weights': self.use_weights.isChecked(),
                'calculate_snr': self.calc_snr.isChecked(),
                'binsize': self.binsize.value(),
            }
            
            # Apply to all items in master table
            item_count = self.controller.master_table.get_item_count()
            if item_count == 0:
                QMessageBox.warning(self, "No Items", "No batch items found to apply settings to.")
                return
            
            # Update all items
            for row_index in range(item_count):
                self.controller.master_table.update_analysis(
                    row_index,
                    continuum_method=settings['continuum_method'],
                    continuum_order=settings['continuum_order'],
                    use_weights=settings['use_weights'],
                    calculate_snr=settings['calculate_snr'],
                    binsize=settings['binsize']
                )
            
            # Update status
            self.settings_status.setText(f"Settings applied to {item_count} items.")
            self.settings_status.setStyleSheet("QLabel { color: green; }")
            
            # Show summary
            method_name = "Polynomial" if continuum_method == "polynomial" else "Flat"
            summary = f"Applied: {method_name} continuum, SNR={settings['calculate_snr']}"
            self.controller.status_updated.emit(summary)
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to apply settings: {str(e)}")
            self.settings_status.setText("Error applying settings.")
            self.settings_status.setStyleSheet("QLabel { color: red; }")
    
    def run_batch(self):
        """Run batch processing on all items."""
        # Check if settings have been applied
        if "not applied" in self.settings_status.text() or "changed" in self.settings_status.text():
            reply = QMessageBox.question(
                self, "Settings Not Applied",
                "Processing settings have not been applied to the batch items. "
                "Do you want to apply current settings before processing?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            )
            
            if reply == QMessageBox.Yes:
                self.apply_settings_to_master_table()
            elif reply == QMessageBox.Cancel:
                return
            # If No, continue with existing settings in master table
        
        self._start_processing()
    
    def _start_processing(self):
        """Start the batch processing thread."""
        
        # Create processor and thread
        self.processor = BatchProcessor(self.controller)
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
        self.apply_settings_btn.setEnabled(False)
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
        self.apply_settings_btn.setEnabled(True)
        
        # Emit signal that batch processing completed
        self.batch_completed.emit(success, results)
    
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
        self.error_log.append(f"Error processing item {index}: {error}")
    
    def cancel_processing(self):
        """Cancel the current batch processing."""
        if self.processor:
            self.processor.cancel()
        if self.processor_thread and self.processor_thread.isRunning():
            self.processor_thread.quit()
            self.processor_thread.wait()
        
        self.run_btn.setEnabled(True)
        self.apply_settings_btn.setEnabled(True)
        self.status_label.setText("Processing cancelled.")