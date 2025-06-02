# rbcodes/GUIs/specgui/batch/panels/processing_panel.py
import os
from datetime import datetime
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
        self.selected_items = []
    
    def process_selected(self, selected_items):
        """Process only the selected batch items."""
        success, results = self.controller.process_batch(selected_items)
        self.finished.emit(success, results)
    
    def process(self):
        """Process batch items - kept for backward compatibility."""
        success, results = self.controller.process_batch(None)
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
        self.selected_from_review = []  # Track items selected from review panel
        
        self.init_ui()
        
        # Connect to controller signals
        self.controller.batch_progress.connect(self.update_progress)
        self.controller.batch_item_completed.connect(self.on_item_completed)
        self.controller.batch_item_failed.connect(self.on_item_failed)
        
        # Connect to table changes to update item counts
        self.controller.table_changed.connect(self.update_process_options)
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Create scroll area
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        
        # Create content widget for scroll area
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        
        # Instructions
        instructions = QLabel(
            "Configure processing settings and choose which items to process. "
            "Click 'Apply Settings' to save your choices, then 'Run Batch Process' to start."
        )
        instructions.setWordWrap(True)
        content_layout.addWidget(instructions)
        
        # Item selection section
        selection_group = QGroupBox("Item Selection")
        selection_layout = QFormLayout()
        
        # Process options dropdown
        self.process_option = QComboBox()
        self.process_option.currentTextChanged.connect(self.on_process_option_changed)
        selection_layout.addRow("Process:", self.process_option)
        
        # Info label showing what will be processed
        self.selection_info = QLabel("Ready to configure item selection.")
        self.selection_info.setWordWrap(True)
        self.selection_info.setStyleSheet("QLabel { color: #666; font-style: italic; }")
        selection_layout.addRow("", self.selection_info)
        
        selection_group.setLayout(selection_layout)
        content_layout.addWidget(selection_group)
        
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
        self.apply_settings_btn = QPushButton("Apply Settings to Selected Items")
        self.apply_settings_btn.clicked.connect(self.apply_settings_to_selected_items)
        self.apply_settings_btn.setToolTip("Save these settings and apply them to the selected items")
        settings_layout.addRow("", self.apply_settings_btn)
        
        # Status for settings
        self.settings_status = QLabel("Settings not applied yet.")
        self.settings_status.setStyleSheet("QLabel { color: red; font-style: italic; }")
        settings_layout.addRow("Status:", self.settings_status)
        
        settings_group.setLayout(settings_layout)
        content_layout.addWidget(settings_group)
        
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
        self.run_btn.setMaximumWidth(200)
        
        # Center the button
        button_container = QHBoxLayout()
        button_container.addStretch()
        button_container.addWidget(self.run_btn)
        button_container.addStretch()
        
        progress_layout.addLayout(button_container)
        
        progress_group.setLayout(progress_layout)
        content_layout.addWidget(progress_group)
        
        # Error log
        error_group = QGroupBox("Error Log")
        error_layout = QVBoxLayout()
        
        self.error_log = QTextEdit()
        self.error_log.setReadOnly(True)
        error_layout.addWidget(self.error_log)
        
        error_group.setLayout(error_layout)
        content_layout.addWidget(error_group)
        
        # Set the content widget to the scroll area
        scroll_area.setWidget(content_widget)
        
        # Add scroll area to main layout
        main_layout.addWidget(scroll_area)
        
        # Initialize selection tracking and update options
        self.update_process_options()
    
    def on_cont_method_changed(self, method):
        """Handle changes to the continuum fitting method."""
        # Enable/disable polynomial order based on method
        is_polynomial = method == "Polynomial"
        self.poly_order.setEnabled(is_polynomial)
        
        # Mark settings as not applied
        self.settings_status.setText("Settings changed - click 'Apply Settings' to save.")
        self.settings_status.setStyleSheet("QLabel { color: red; font-style: italic; }")
    
    def update_process_options(self):
        """Update the process options dropdown with current item counts."""
        # Get current selection to preserve it
        current_selection = self.process_option.currentText()
        
        # Clear and repopulate options
        self.process_option.clear()
        
        # Get item counts by status
        ready_count = len(self.controller.master_table.get_items_by_status("ready"))
        complete_count = len(self.controller.master_table.get_items_by_status("complete"))
        error_count = len(self.controller.master_table.get_items_by_status("error"))
        needs_count = len(self.controller.master_table.get_items_by_status("needs_processing"))
        total_count = self.controller.get_item_count()
        
        # Build option list with counts
        options = []
        
        # Ready items (default)
        incomplete_count = ready_count + error_count + needs_count
        if incomplete_count > 0:
            options.append(f"Ready Items Only ({incomplete_count} items)")
        
        # All items
        if total_count > 0:
            options.append(f"All Items ({total_count} items)")
        
        # Failed items only
        if error_count > 0:
            options.append(f"Failed Items Only ({error_count} items)")
        
        # Selected from review
        if self.selected_from_review:
            options.append(f"Selected from Review ({len(self.selected_from_review)} items)")
        
        # Add options to dropdown
        if options:
            self.process_option.addItems(options)
            
            # Try to restore previous selection
            for i, option in enumerate(options):
                if current_selection and current_selection.split("(")[0].strip() in option:
                    self.process_option.setCurrentIndex(i)
                    break
            else:
                # Default to first option (Ready Items Only)
                self.process_option.setCurrentIndex(0)
        else:
            self.process_option.addItem("No items to process")
        
        # Update info label
        self.on_process_option_changed(self.process_option.currentText())
    
    def on_process_option_changed(self, option_text):
        """Handle changes to the process option selection."""
        if "No items" in option_text:
            self.selection_info.setText("Add batch items in the Configuration tab first.")
            self.run_btn.setEnabled(False)
            return
        
        # Extract the selection type and count
        if "Ready Items Only" in option_text:
            self.selection_info.setText("Will process items with status: ready, error, or needs_processing.")
        elif "All Items" in option_text:
            self.selection_info.setText("Will process ALL items, including already completed ones.")
        elif "Failed Items Only" in option_text:
            self.selection_info.setText("Will process only items that failed in previous runs.")
        elif "Selected from Review" in option_text:
            self.selection_info.setText("Will process only the items you selected in the Review panel.")
        
        self.run_btn.setEnabled(True)
        
        # Mark settings as needing reapplication if selection changed
        if "not applied" not in self.settings_status.text():
            self.settings_status.setText("Selection changed - click 'Apply Settings' to update.")
            self.settings_status.setStyleSheet("QLabel { color: orange; font-style: italic; }")
    
    def set_selected_from_review(self, selected_indices):
        """Set the items selected from the review panel."""
        self.selected_from_review = selected_indices
        self.update_process_options()
        
        # Auto-select the "Selected from Review" option if items were provided
        if selected_indices:
            for i in range(self.process_option.count()):
                if "Selected from Review" in self.process_option.itemText(i):
                    self.process_option.setCurrentIndex(i)
                    break
    
    def get_items_to_process(self):
        """Get the list of row indices to process based on current selection."""
        option_text = self.process_option.currentText()
        
        if "Ready Items Only" in option_text:
            # Get items that need processing
            ready_items = self.controller.master_table.get_items_by_status("ready")
            error_items = self.controller.master_table.get_items_by_status("error")
            needs_items = self.controller.master_table.get_items_by_status("needs_processing")
            return ready_items + error_items + needs_items
        
        elif "All Items" in option_text:
            # Process everything
            return list(range(self.controller.get_item_count()))
        
        elif "Failed Items Only" in option_text:
            # Only failed items
            return self.controller.master_table.get_items_by_status("error")
        
        elif "Selected from Review" in option_text:
            # Use the stored selection from review panel
            return self.selected_from_review.copy()
        
        else:
            # Fallback - no items
            return []
    
    def apply_settings_to_selected_items(self):
        """Apply current UI settings to the selected items only."""
        try:
            # Get the items that will be processed
            items_to_process = self.get_items_to_process()
            
            if not items_to_process:
                QMessageBox.warning(self, "No Items", "No items selected for processing.")
                return
            
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
            
            # Apply to selected items only
            for row_index in items_to_process:
                self.controller.master_table.update_analysis(
                    row_index,
                    continuum_method=settings['continuum_method'],
                    continuum_order=settings['continuum_order'],
                    use_weights=settings['use_weights'],
                    calculate_snr=settings['calculate_snr'],
                    binsize=settings['binsize']
                )
            
            # Update status
            option_text = self.process_option.currentText()
            self.settings_status.setText(f"Settings applied to {len(items_to_process)} items ({option_text.split('(')[0].strip()}).")
            self.settings_status.setStyleSheet("QLabel { color: green; }")
            
            # Show summary
            method_name = "Polynomial" if continuum_method == "polynomial" else "Flat"
            summary = f"Applied {method_name} continuum settings to {len(items_to_process)} selected items"
            self.controller.status_updated.emit(summary)
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to apply settings: {str(e)}")
            self.settings_status.setText("Error applying settings.")
            self.settings_status.setStyleSheet("QLabel { color: red; }")
    
    def run_batch(self):
        """Run batch processing on the selected items."""
        # Check if settings have been applied
        if "not applied" in self.settings_status.text() or "changed" in self.settings_status.text():
            reply = QMessageBox.question(
                self, "Settings Not Applied",
                "Processing settings have not been applied to the selected items. "
                "Do you want to apply current settings before processing?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            )
            
            if reply == QMessageBox.Yes:
                self.apply_settings_to_selected_items()
            elif reply == QMessageBox.Cancel:
                return
            # If No, continue with existing settings in master table
        
        # Get items to process
        items_to_process = self.get_items_to_process()
        
        if not items_to_process:
            QMessageBox.warning(self, "No Items", "No items selected for processing.")
            return
        
        # Show confirmation for what will be processed
        option_text = self.process_option.currentText()
        reply = QMessageBox.question(
            self, "Confirm Processing",
            f"Start batch processing?\n\nWill process: {option_text}\n\nThis may take several minutes.",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.No:
            return
        
        self._start_processing(items_to_process)
    
    def _start_processing(self, selected_items):
        """Start the batch processing thread with selected items."""
        
        # Create processor and thread
        self.processor = BatchProcessor(self.controller)
        self.processor_thread = QThread()
        self.processor.moveToThread(self.processor_thread)
        
        # Set the selected items for processing
        self.processor.selected_items = selected_items
        
        # Connect signals
        self.processor_thread.started.connect(lambda: self.processor.process_selected(selected_items))
        self.processor.progress.connect(self.update_progress)
        self.processor.item_completed.connect(self.on_item_completed)
        self.processor.item_failed.connect(self.on_item_failed)
        self.processor.finished.connect(self._on_processing_finished)
        
        # Update UI
        self.run_btn.setEnabled(False)
        self.apply_settings_btn.setEnabled(False)
        self.progress_bar.setValue(0)
        self.error_log.clear()
        
        # Show what's being processed
        option_text = self.process_option.currentText()
        self.controller.status_updated.emit(f"Starting batch processing: {option_text}")
        
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