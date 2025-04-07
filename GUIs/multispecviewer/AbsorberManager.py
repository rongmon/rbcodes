import pandas as pd
from PyQt5.QtWidgets import (QWidget, QTableWidget, QTableWidgetItem, 
                             QHBoxLayout, QVBoxLayout, QComboBox, QPushButton,
                             QHeaderView, QSizePolicy)
from PyQt5.QtGui import QColor
from PyQt5 import QtCore, QtGui
from utils import rb_utility as rt
clr = rt.rb_set_color()

class AbsorberManager(QWidget):
    """
    A widget for managing multiple absorber systems with their redshifts, line lists, and colors.
    Allows adding, removing, plotting, and hiding absorbers.
    """
    
    def __init__(self, parent=None, colors=None):
        super().__init__(parent)
        # Set up the dark theme after creating the widgets
        self.setup_dark_theme()
        self.parent = parent
        
        # Import color utility and use those colors if not provided
        self.colors = colors or list(clr.keys())[1:]  # Skip the first color (usually background)        

        # Define available line lists
        self.line_options = ['LLS', 'LLS Small', 'DLA', 'LBG', 'Gal', 'Eiger_Strong', 'None']
        
        # Initialize data storage for absorbers
        self.absorbers_df = pd.DataFrame(columns=['Zabs', 'LineList', 'Color'])
        
        # Initialize UI components
        self.initUI()
        
    def initUI(self):
        """Initialize the user interface components"""
        self.layout = QVBoxLayout(self)
        
        # Create table for displaying absorbers
        self.table = QTableWidget(1, 6)  # Start with 1 row, will expand as needed
        self.table.setHorizontalHeaderLabels(['Line Lists', 'z', 'Color', 'Plot', 'Remove', 'Hide'])
        
        # Configure table properties
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Set default widgets for initial rows
        self._add_row_widgets(0)
        self._add_row_widgets(1)
        
        # Add table to layout
        self.layout.addWidget(self.table)
        
        # Connect cell changed signal
        self.table.cellChanged.connect(self.on_cell_changed)
        
    def _add_row_widgets(self, row_index):
        """Add widgets to a specific row in the table"""
        # Line list combo box
        line_combo = QComboBox()
        for item in self.line_options:
            line_combo.addItem(item)
        self.table.setCellWidget(row_index, 0, line_combo)
        
        # Color combo box
        color_combo = QComboBox()
        for color in self.colors:
            color_combo.addItem(color)
        self.table.setCellWidget(row_index, 2, color_combo)
        
        # Plot button
        plot_btn = QPushButton("Plot")
        plot_btn.clicked.connect(lambda: self.plot_absorber(row_index))
        self.table.setCellWidget(row_index, 3, plot_btn)
        
        # Remove button
        remove_btn = QPushButton("Remove")
        remove_btn.clicked.connect(lambda: self.remove_absorber(row_index))
        self.table.setCellWidget(row_index, 4, remove_btn)
        
        # Hide button
        hide_btn = QPushButton("Hide")
        hide_btn.clicked.connect(lambda: self.toggle_hide_absorber(row_index))
        self.table.setCellWidget(row_index, 5, hide_btn)
    
    def add_absorber(self, z_abs, line_list=None, color=None):
        """Add a new absorber to the manager"""
        # Default values if not provided
        line_list = line_list or 'LLS'
        color = color or 'white'
        
        # Add to dataframe
        new_row = pd.Series({'Zabs': z_abs, 'LineList': line_list, 'Color': color})
        self.absorbers_df = self.absorbers_df.append(new_row, ignore_index=True)
        
        # Find the first empty row or add a new one
        row_found = False
        for row in range(self.table.rowCount()):
            if self.table.item(row, 1) is None:
                # This is an empty row, use it
                self._populate_row(row, z_abs, line_list, color)
                row_found = True
                break
        
        if not row_found:
            # No empty rows, add a new one
            new_row_idx = self.table.rowCount()
            self.table.setRowCount(new_row_idx + 1)
            self._add_row_widgets(new_row_idx)
            self._populate_row(new_row_idx, z_abs, line_list, color)
        
        # Return success status
        return True
    
    def _populate_row(self, row, z_abs, line_list, color):
        """Populate a row with absorber data"""
        # Set redshift value
        self.table.setItem(row, 1, QTableWidgetItem(f"{z_abs:.6f}"))
        
        # Set line list selection
        line_combo = self.table.cellWidget(row, 0)
        line_index = line_combo.findText(line_list)
        if line_index >= 0:
            line_combo.setCurrentIndex(line_index)
        
        # Set color selection
        color_combo = self.table.cellWidget(row, 2)
        color_index = color_combo.findText(color)
        if color_index >= 0:
            color_combo.setCurrentIndex(color_index)
    
    def plot_absorber(self, row):
        """Plot the absorber at the specified row"""
        # Check if row exists and has data
        if row >= self.table.rowCount() or self.table.item(row, 1) is None:
            print(f"Error: No absorber data at row {row}")
            return False
        
        try:
            # Get absorber data
            z_abs = float(self.table.item(row, 1).text())
            line_list = self.table.cellWidget(row, 0).currentText()
            color = self.table.cellWidget(row, 2).currentText()
            
            # Call the parent's plotting method if available
            if hasattr(self.parent, 'plot_absorber_lines'):
                success = self.parent.plot_absorber_lines(row, z_abs, line_list, color)
                return success
            else:
                print("Error: Parent does not have plot_absorber_lines method")
                return False
        except Exception as e:
            print(f"Error plotting absorber: {e}")
            return False
    
    def remove_absorber(self, row):
        """Remove the absorber at the specified row"""
        try:
            # Remove from dataframe if exists
            if row < len(self.absorbers_df):
                self.absorbers_df = self.absorbers_df.drop(row).reset_index(drop=True)
            
            # Clear the row in the table
            self.table.removeRow(row)
            
            # Add a new empty row at the end to maintain minimum rows
            new_row_idx = self.table.rowCount()
            self.table.insertRow(new_row_idx)
            self._add_row_widgets(new_row_idx)
            
            # Call parent's method to remove absorber lines if available
            if hasattr(self.parent, 'remove_absorber_lines'):
                self.parent.remove_absorber_lines(row)
                
            return True
        except Exception as e:
            print(f"Error removing absorber: {e}")
            return False
    
    def toggle_hide_absorber(self, row):
        """Toggle visibility of the absorber at the specified row"""
        try:
            # Call parent's method to toggle visibility if available
            if hasattr(self.parent, 'toggle_absorber_visibility'):
                success = self.parent.toggle_absorber_visibility(row)
                
                # Update the Hide button text if successful
                if success:
                    hide_btn = self.table.cellWidget(row, 5)
                    current_text = hide_btn.text()
                    hide_btn.setText("Show" if current_text == "Hide" else "Hide")
                
                return success
            else:
                print("Error: Parent does not have toggle_absorber_visibility method")
                return False
        except Exception as e:
            print(f"Error toggling absorber visibility: {e}")
            return False
    
    def on_cell_changed(self, row, column):
        """Handle manual edits to cells"""
        if column == 1:  # Redshift column
            try:
                z_abs = float(self.table.item(row, 1).text())
                
                # Update dataframe
                if row < len(self.absorbers_df):
                    self.absorbers_df.at[row, 'Zabs'] = z_abs
                
                # Notify parent if needed
                if hasattr(self.parent, 'update_absorber_redshift'):
                    self.parent.update_absorber_redshift(row, z_abs)
            except ValueError:
                print("Invalid redshift value")
                # Revert to previous value if in dataframe
                if row < len(self.absorbers_df):
                    prev_z = self.absorbers_df.at[row, 'Zabs']
                    self.table.item(row, 1).setText(f"{prev_z:.6f}")
    
    def get_absorber_count(self):
        """Return the number of defined absorbers"""
        return len(self.absorbers_df)
    
    def get_absorber_data(self, row):
        """Get the data for a specific absorber"""
        if row < len(self.absorbers_df):
            return self.absorbers_df.iloc[row].to_dict()
        return None

    def setup_dark_theme(self):
        
        """
        Apply a dark theme similar to the main application's color scheme
        """
        # Dark background colors
        dark_background = QColor(53, 53, 53)
        darker_background = QColor(25, 25, 25)
        
        # Text and accent colors
        white_text = QColor(255, 255, 255)
        highlight_color = QColor(42, 130, 218)
        
        # Fix around line 252-300
        dark_stylesheet = """
        QWidget {
            background-color: """ + dark_background.name() + """;
            color: """ + white_text.name() + """;
            selection-background-color: """ + highlight_color.name() + """;
        }
        
        QTableWidget {
            background-color: """ + darker_background.name() + """;
            alternate-background-color: """ + dark_background.name() + """;
            color: """ + white_text.name() + """;
            selection-background-color: """ + highlight_color.name() + """;
        }
        
        QTableWidget::item {
            background-color: """ + darker_background.name() + """;
            color: """ + white_text.name() + """;
        }
        
        QTableWidget::item:selected {
            background-color: """ + highlight_color.name() + """;
            color: """ + white_text.name() + """;
        }
        
        QHeaderView::section {
            background-color: """ + dark_background.name() + """;
            color: """ + white_text.name() + """;
            padding: 5px;
            border: 1px solid """ + darker_background.name() + """;
        }
        
        QPushButton {
            background-color: """ + dark_background.name() + """;
            color: """ + white_text.name() + """;
            border: 1px solid """ + highlight_color.name() + """;
            padding: 5px;
        }
        
        QPushButton:hover {
            background-color: """ + highlight_color.name() + """;
            color: """ + white_text.name() + """;
        }
        
        QComboBox {
            background-color: """ + darker_background.name() + """;
            color: """ + white_text.name() + """;
            selection-background-color: """ + highlight_color.name() + """;
        }
        
        QComboBox::drop-down {
            background-color: """ + dark_background.name() + """;
        }
        """
        # Apply stylesheet
        self.setStyleSheet(dark_stylesheet)
        
        # Try to set alternating row colors for all table widgets
        for child in self.findChildren(QTableWidget):
            child.setAlternatingRowColors(True)
            child.setStyleSheet("""
                QTableWidget { 
                    alternate-background-color: #353535; 
                    background-color: #191919; 
                }
            """)