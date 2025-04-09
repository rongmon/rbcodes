#!/usr/bin/env python
"""
Absorber Manager Widget - UI component for managing absorption systems.
"""
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem,
    QPushButton, QComboBox, QHeaderView
)
from PyQt5.QtCore import pyqtSignal

from rbcodes.GUIs.spectral_gui.utils.constants import COLORS, LINE_LISTS


class AbsorberManagerWidget(QWidget):
    """
    Widget for managing absorption systems.
    
    This widget provides a table interface for viewing, editing, and
    interacting with absorbers.
    
    Parameters
    ----------
    absorber_manager : AbsorberManager
        Manager containing absorber data
    parent : QWidget, optional
        Parent widget, by default None
    """
    
    # Signals for absorber actions
    plotRequested = pyqtSignal(int)
    hideRequested = pyqtSignal(int)
    removeRequested = pyqtSignal(int)
    
    def __init__(self, absorber_manager, parent=None):
        """Initialize the absorber manager widget."""
        super(AbsorberManagerWidget, self).__init__(parent)
        
        # Store reference to absorber manager
        self.absorber_manager = absorber_manager
        
        # Set up the UI
        self._create_ui()
        
        # Fill table with current data
        self.refresh()
    
    def _create_ui(self):
        """Create the user interface."""
        # Create layout
        layout = QVBoxLayout(self)
        
        # Create table
        self.table = QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(
            ['Line Lists', 'z', 'Color', 'Plot', 'Remove', 'Hide']
        )
        
        # Configure table
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeToContents)
        
        # Connect cell changed signal
        self.table.cellChanged.connect(self._on_cell_changed)
        
        # Add table to layout
        layout.addWidget(self.table)
        
        # Set layout
        self.setLayout(layout)
    
    def refresh(self):
        """Refresh the table with current data."""
        # Disconnect cell changed signal temporarily
        self.table.cellChanged.disconnect(self._on_cell_changed)
        
        # Clear existing rows
        self.table.setRowCount(0)
        
        # Add rows for each absorber
        for i, absorber in enumerate(self.absorber_manager.absorbers):
            self._add_absorber_row(i, absorber)
        
        # Add an empty row at the end
        self._add_empty_row()
        
        # Reconnect cell changed signal
        self.table.cellChanged.connect(self._on_cell_changed)
    
    def _add_absorber_row(self, index, absorber):
        """
        Add a row for an absorber.
        
        Parameters
        ----------
        index : int
            Index of the absorber
        absorber : Absorber
            Absorber object
        """
        # Add a new row
        row = self.table.rowCount()
        self.table.insertRow(row)
        
        # Line list combo
        line_list_combo = QComboBox()
        for line_list in LINE_LISTS:
            line_list_combo.addItem(line_list)
        line_list_combo.setCurrentText(absorber.line_list)
        line_list_combo.currentIndexChanged.connect(
            lambda: self._on_line_list_changed(row)
        )
        self.table.setCellWidget(row, 0, line_list_combo)
        
        # Redshift
        redshift_item = QTableWidgetItem(f"{absorber.redshift:.4f}")
        self.table.setItem(row, 1, redshift_item)
        
        # Color combo
        color_combo = QComboBox()
        for color in COLORS:
            color_combo.addItem(color)
        color_combo.setCurrentText(absorber.color)
        color_combo.currentIndexChanged.connect(
            lambda: self._on_color_changed(row)
        )
        self.table.setCellWidget(row, 2, color_combo)
        
        # Plot button
        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(lambda: self._on_plot_clicked(row))
        self.table.setCellWidget(row, 3, plot_button)
        
        # Remove button
        remove_button = QPushButton("Remove")
        remove_button.clicked.connect(lambda: self._on_remove_clicked(row))
        self.table.setCellWidget(row, 4, remove_button)
        
        # Hide button
        hide_button = QPushButton("Hide")
        hide_button.clicked.connect(lambda: self._on_hide_clicked(row))
        self.table.setCellWidget(row, 5, hide_button)
    
    def _add_empty_row(self):
        """Add an empty row at the end of the table."""
        # Add a new row
        row = self.table.rowCount()
        self.table.insertRow(row)
        
        # Line list combo
        line_list_combo = QComboBox()
        for line_list in LINE_LISTS:
            line_list_combo.addItem(line_list)
        self.table.setCellWidget(row, 0, line_list_combo)
        
        # Color combo
        color_combo = QComboBox()
        for color in COLORS:
            color_combo.addItem(color)
        self.table.setCellWidget(row, 2, color_combo)
        
        # Plot button
        plot_button = QPushButton("Plot")
        plot_button.setEnabled(False)
        self.table.setCellWidget(row, 3, plot_button)
        
        # Remove button
        remove_button = QPushButton("Remove")
        remove_button.setEnabled(False)
        self.table.setCellWidget(row, 4, remove_button)
        
        # Hide button
        hide_button = QPushButton("Hide")
        hide_button.setEnabled(False)
        self.table.setCellWidget(row, 5, hide_button)
    
    def _on_cell_changed(self, row, column):
        """
        Handle cell value changes.
        
        Parameters
        ----------
        row : int
            Row index
        column : int
            Column index
        """
        # Only handle redshift changes (column 1)
        if column != 1:
            return
        
        # Get the new value
        item = self.table.item(row, column)
        if not item:
            return
        
        # Try to parse as float
        try:
            redshift = float(item.text())
            
            # Update the absorber if row is valid
            if row < len(self.absorber_manager.absorbers):
                absorber = self.absorber_manager.get_absorber(row)
                if absorber:
                    absorber.redshift = redshift
                            # Otherwise create a new absorber
            elif item.text().strip():
                from rbcodes.GUIs.spectral_gui.core.absorber import Absorber
                
                # Get line list and color
                line_list_combo = self.table.cellWidget(row, 0)
                color_combo = self.table.cellWidget(row, 2)
                
                # Create absorber
                absorber = Absorber(
                    redshift=redshift,
                    line_list=line_list_combo.currentText(),
                    color=color_combo.currentText()
                )
                
                # Add to manager
                self.absorber_manager.add_absorber(absorber)
                
                # Enable buttons
                self.table.cellWidget(row, 3).setEnabled(True)
                self.table.cellWidget(row, 4).setEnabled(True)
                self.table.cellWidget(row, 5).setEnabled(True)
                
                # Add a new empty row
                self._add_empty_row()
        except ValueError:
            # Reset to original value if invalid
            if row < len(self.absorber_manager.absorbers):
                absorber = self.absorber_manager.get_absorber(row)
                if absorber:
                    item.setText(f"{absorber.redshift:.4f}")
            else:
                item.setText("")
    
    def _on_line_list_changed(self, row):
        """
        Handle line list selection change.
        
        Parameters
        ----------
        row : int
            Row index
        """
        # Check if row is valid
        if row >= len(self.absorber_manager.absorbers):
            return
        
        # Get the selected line list
        line_list_combo = self.table.cellWidget(row, 0)
        line_list = line_list_combo.currentText()
        
        # Update the absorber
        absorber = self.absorber_manager.get_absorber(row)
        if absorber:
            absorber.line_list = line_list
    
    def _on_color_changed(self, row):
        """
        Handle color selection change.
        
        Parameters
        ----------
        row : int
            Row index
        """
        # Check if row is valid
        if row >= len(self.absorber_manager.absorbers):
            return
        
        # Get the selected color
        color_combo = self.table.cellWidget(row, 2)
        color = color_combo.currentText()
        
        # Update the absorber
        absorber = self.absorber_manager.get_absorber(row)
        if absorber:
            absorber.color = color
    
    def _on_plot_clicked(self, row):
        """
        Handle plot button click.
        
        Parameters
        ----------
        row : int
            Row index
        """
        # Check if row is valid
        if row >= len(self.absorber_manager.absorbers):
            return
        
        # Emit signal
        self.plotRequested.emit(row)
    
    def _on_remove_clicked(self, row):
        """
        Handle remove button click.
        
        Parameters
        ----------
        row : int
            Row index
        """
        # Check if row is valid
        if row >= len(self.absorber_manager.absorbers):
            return
        
        # Emit signal
        self.removeRequested.emit(row)
    
    def _on_hide_clicked(self, row):
        """
        Handle hide button click.
        
        Parameters
        ----------
        row : int
            Row index
        """
        # Check if row is valid
        if row >= len(self.absorber_manager.absorbers):
            return
        
        # Get the button
        button = self.table.cellWidget(row, 5)
        
        # Toggle button state
        if button.text() == "Hide":
            button.setText("Show")
            button.setStyleSheet("background-color: green")
            self.hideRequested.emit(row)
        else:
            button.setText("Hide")
            button.setStyleSheet("")
            self.plotRequested.emit(row)