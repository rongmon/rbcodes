# rbcodes/GUIs/specgui/batch/panels/batch_selection_dialog.py
"""
Batch Selection Dialog for choosing multiple systems for processing.

This dialog provides a familiar scrollable checklist interface with smart
selection buttons for choosing systems to process in batch mode.
"""

import os
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QTableWidget, QTableWidgetItem,
                           QHeaderView, QCheckBox, QGroupBox, QDialogButtonBox,
                           QMessageBox, QLineEdit, QSplitter)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor


class BatchSelectionDialog(QDialog):
    """Dialog for selecting multiple systems for batch processing."""
    
    def __init__(self, controller, parent=None):
        super().__init__(parent)
        self.controller = controller
        self.selected_indices = []
        
        self.setWindowTitle("Select Systems for Processing")
        self.setMinimumSize(800, 600)
        self.setModal(True)
        
        self.init_ui()
        self.populate_systems()
    
    def init_ui(self):
        """Initialize the user interface."""
        layout = QVBoxLayout(self)
        
        # Instructions
        instructions = QLabel(
            "Select systems to process using the buttons below or by checking individual items. "
            "You can combine smart selection with manual adjustments."
        )
        instructions.setWordWrap(True)
        layout.addWidget(instructions)
        
        # Smart selection buttons
        selection_group = QGroupBox("Quick Selection")
        selection_layout = QVBoxLayout()
        
        # Row 1: Status-based selection
        status_buttons = QHBoxLayout()
        
        self.select_failed_btn = QPushButton("Select All Failed")
        self.select_failed_btn.clicked.connect(self.select_failed_items)
        self.select_failed_btn.setToolTip("Select all systems with error status")
        
        self.select_pending_btn = QPushButton("Select All Pending")
        self.select_pending_btn.clicked.connect(self.select_pending_items)
        self.select_pending_btn.setToolTip("Select all systems ready for processing")
        
        self.select_complete_btn = QPushButton("Select All Complete")
        self.select_complete_btn.clicked.connect(self.select_complete_items)
        self.select_complete_btn.setToolTip("Select all successfully processed systems")
        
        status_buttons.addWidget(self.select_failed_btn)
        status_buttons.addWidget(self.select_pending_btn)
        status_buttons.addWidget(self.select_complete_btn)
        status_buttons.addStretch()
        
        # Row 2: General selection
        general_buttons = QHBoxLayout()
        
        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.clicked.connect(self.select_all_items)
        
        self.select_none_btn = QPushButton("Clear Selection")
        self.select_none_btn.clicked.connect(self.clear_selection)
        
        self.invert_selection_btn = QPushButton("Invert Selection")
        self.invert_selection_btn.clicked.connect(self.invert_selection)
        
        general_buttons.addWidget(self.select_all_btn)
        general_buttons.addWidget(self.select_none_btn)
        general_buttons.addWidget(self.invert_selection_btn)
        general_buttons.addStretch()
        
        selection_layout.addLayout(status_buttons)
        selection_layout.addLayout(general_buttons)
        selection_group.setLayout(selection_layout)
        layout.addWidget(selection_group)
        
        # Systems table
        table_group = QGroupBox("Systems")
        table_layout = QVBoxLayout()
        
        # Selection counter
        self.selection_counter = QLabel("0 systems selected for processing")
        self.selection_counter.setStyleSheet("QLabel { font-weight: bold; color: #0066cc; }")
        table_layout.addWidget(self.selection_counter)
        
        # Table
        self.systems_table = QTableWidget()
        self.systems_table.setColumnCount(8)
        self.systems_table.setHorizontalHeaderLabels([
            "", "System", "Filename", "Transition", "Redshift", "EW Range", "Status", "Results"
        ])
        
        # Set column widths
        header = self.systems_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Fixed)  # Checkbox
        header.resizeSection(0, 30)
        header.setSectionResizeMode(1, QHeaderView.Fixed)  # System number
        header.resizeSection(1, 60)
        header.setSectionResizeMode(2, QHeaderView.Stretch)  # Filename
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)  # Transition
        header.setSectionResizeMode(4, QHeaderView.ResizeToContents)  # Redshift
        header.setSectionResizeMode(5, QHeaderView.ResizeToContents)  # EW Range
        header.setSectionResizeMode(6, QHeaderView.ResizeToContents)  # Status
        header.setSectionResizeMode(7, QHeaderView.ResizeToContents)  # Results
        
        self.systems_table.setRowCount(0)
        self.systems_table.setAlternatingRowColors(True)
        self.systems_table.setSelectionBehavior(QTableWidget.SelectRows)
        
        # Connect checkbox changes to update counter
        self.systems_table.itemChanged.connect(self.update_selection_counter)
        
        table_layout.addWidget(self.systems_table)
        table_group.setLayout(table_layout)
        layout.addWidget(table_group)
        
        # Dialog buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept_selection)
        button_box.rejected.connect(self.reject)
        
        # Add preview button
        self.preview_btn = QPushButton("Preview Selection")
        self.preview_btn.clicked.connect(self.preview_selection)
        button_box.addButton(self.preview_btn, QDialogButtonBox.ActionRole)
        
        layout.addWidget(button_box)
    
    def populate_systems(self):
        """Populate the systems table from master table."""
        self.systems_table.setRowCount(0)
        
        # Get all items from master table
        items = self.controller.master_table.get_all_items()
        
        if not items:
            return
        
        # Populate table
        for i, item in enumerate(items):
            row = self.systems_table.rowCount()
            self.systems_table.insertRow(row)
            
            # Checkbox
            checkbox_item = QTableWidgetItem()
            checkbox_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            checkbox_item.setCheckState(Qt.Unchecked)
            self.systems_table.setItem(row, 0, checkbox_item)
            
            # System number
            self.systems_table.setItem(row, 1, QTableWidgetItem(f"{i+1}"))
            
            # Filename
            filename = os.path.basename(item.template.filename)
            self.systems_table.setItem(row, 2, QTableWidgetItem(filename))
            
            # Transition
            transition_name = item.template.transition_name
            transition_wave = item.template.transition
            transition_text = f"{transition_name} ({transition_wave:.2f} Å)"
            self.systems_table.setItem(row, 3, QTableWidgetItem(transition_text))
            
            # Redshift
            self.systems_table.setItem(row, 4, QTableWidgetItem(f"{item.template.redshift:.6f}"))
            
            # EW Range
            ew_range_text = f"[{item.template.ew_vmin}, {item.template.ew_vmax}]"
            self.systems_table.setItem(row, 5, QTableWidgetItem(ew_range_text))
            
            # Status with color
            status = item.analysis.processing_status
            status_text = status.replace('_', ' ').title()
            status_item = QTableWidgetItem(status_text)
            
            if status == 'complete':
                status_item.setBackground(QColor(200, 255, 200))
            elif status == 'error':
                status_item.setBackground(QColor(255, 200, 200))
            elif 'processing' in status:
                status_item.setBackground(QColor(255, 255, 200))
            elif status in ['needs_processing', 'needs_ew_recalc', 'ready']:
                status_item.setBackground(QColor(255, 240, 200))
            
            self.systems_table.setItem(row, 6, status_item)
            
            # Results
            if hasattr(item.results, 'W') and item.results.W > 0:
                ew_value = item.results.W
                ew_error = item.results.W_e
                results_text = f"EW: {ew_value:.3f}±{ew_error:.3f}"
            else:
                results_text = "Not calculated"
            
            self.systems_table.setItem(row, 7, QTableWidgetItem(results_text))
        
        # Update counter
        self.update_selection_counter()
    
    def update_selection_counter(self):
        """Update the selection counter label."""
        selected_count = 0
        total_count = self.systems_table.rowCount()
        
        for row in range(total_count):
            checkbox_item = self.systems_table.item(row, 0)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                selected_count += 1
        
        self.selection_counter.setText(f"{selected_count} of {total_count} systems selected for processing")
        
        # Update button text based on selection
        ok_button = self.findChild(QDialogButtonBox).button(QDialogButtonBox.Ok)
        if ok_button:
            if selected_count > 0:
                ok_button.setText(f"Process {selected_count} Systems")
            else:
                ok_button.setText("OK")
    
    def select_failed_items(self):
        """Select all items with error status."""
        items = self.controller.master_table.get_all_items()
        
        for row, item in enumerate(items):
            checkbox_item = self.systems_table.item(row, 0)
            if checkbox_item:
                if item.analysis.processing_status == 'error':
                    checkbox_item.setCheckState(Qt.Checked)
        
        self.update_selection_counter()
    
    def select_pending_items(self):
        """Select all items ready for processing."""
        items = self.controller.master_table.get_all_items()
        
        for row, item in enumerate(items):
            checkbox_item = self.systems_table.item(row, 0)
            if checkbox_item:
                status = item.analysis.processing_status
                if status in ['ready', 'needs_processing', 'needs_ew_recalc']:
                    checkbox_item.setCheckState(Qt.Checked)
        
        self.update_selection_counter()
    
    def select_complete_items(self):
        """Select all completed items."""
        items = self.controller.master_table.get_all_items()
        
        for row, item in enumerate(items):
            checkbox_item = self.systems_table.item(row, 0)
            if checkbox_item:
                if item.analysis.processing_status == 'complete':
                    checkbox_item.setCheckState(Qt.Checked)
        
        self.update_selection_counter()
    
    def select_all_items(self):
        """Select all items."""
        for row in range(self.systems_table.rowCount()):
            checkbox_item = self.systems_table.item(row, 0)
            if checkbox_item:
                checkbox_item.setCheckState(Qt.Checked)
        
        self.update_selection_counter()
    
    def clear_selection(self):
        """Clear all selections."""
        for row in range(self.systems_table.rowCount()):
            checkbox_item = self.systems_table.item(row, 0)
            if checkbox_item:
                checkbox_item.setCheckState(Qt.Unchecked)
        
        self.update_selection_counter()
    
    def invert_selection(self):
        """Invert the current selection."""
        for row in range(self.systems_table.rowCount()):
            checkbox_item = self.systems_table.item(row, 0)
            if checkbox_item:
                current_state = checkbox_item.checkState()
                new_state = Qt.Checked if current_state == Qt.Unchecked else Qt.Unchecked
                checkbox_item.setCheckState(new_state)
        
        self.update_selection_counter()
    
    def preview_selection(self):
        """Show a preview of the selected systems."""
        selected_indices = self.get_selected_indices()
        
        if not selected_indices:
            QMessageBox.information(self, "Preview", "No systems selected.")
            return
        
        # Build preview text
        preview_text = f"Selected {len(selected_indices)} systems for processing:\n\n"
        
        items = self.controller.master_table.get_all_items()
        for i, index in enumerate(selected_indices):
            if index < len(items):
                item = items[index]
                filename = os.path.basename(item.template.filename)
                transition = item.template.transition_name
                status = item.analysis.processing_status.replace('_', ' ').title()
                preview_text += f"{i+1}. System {index+1}: {filename} | {transition} | {status}\n"
        
        # Show preview dialog
        QMessageBox.information(self, "Selection Preview", preview_text)
    
    def get_selected_indices(self):
        """Get the indices of selected systems."""
        selected_indices = []
        
        for row in range(self.systems_table.rowCount()):
            checkbox_item = self.systems_table.item(row, 0)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                selected_indices.append(row)
        
        return selected_indices
    
    def accept_selection(self):
        """Accept the current selection."""
        selected_indices = self.get_selected_indices()
        
        if not selected_indices:
            reply = QMessageBox.question(
                self, "No Selection", 
                "No systems are selected. Do you want to close without processing anything?",
                QMessageBox.Yes | QMessageBox.No
            )
            
            if reply == QMessageBox.No:
                return
        
        self.selected_indices = selected_indices
        self.accept()
    
    def get_selection(self):
        """Get the final selection after dialog closes."""
        return self.selected_indices


def show_batch_selection_dialog(controller, parent=None):
    """
    Convenience function to show the batch selection dialog.
    
    Returns:
        list: Selected system indices, or empty list if cancelled
    """
    dialog = BatchSelectionDialog(controller, parent)
    
    if dialog.exec_() == QDialog.Accepted:
        return dialog.get_selection()
    else:
        return []