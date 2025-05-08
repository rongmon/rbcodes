import pandas as pd
from PyQt5.QtWidgets import (QWidget, QTableWidget, QTableWidgetItem, 
                             QHBoxLayout, QVBoxLayout, QComboBox, QPushButton,
                             QHeaderView, QSizePolicy, QApplication, QCheckBox)
from PyQt5.QtGui import QColor
from PyQt5 import QtCore, QtGui, QtWidgets
from rbcodes.utils import rb_utility as rt
clr = rt.rb_set_color()

class AbsorberManager(QWidget):
    """
    A widget for managing multiple absorber systems with their redshifts, line lists, and colors.
    Allows adding, removing, plotting, and hiding absorbers using a unified checkbox approach.
    
    Note: To fully support this widget, the parent should implement:
    - plot_absorber_lines(row, z_abs, line_list, color) - Plot lines for this absorber
    - remove_absorber_lines(row) - Remove lines for this absorber
    - toggle_absorber_visibility(row, visible) - Toggle or set visibility of absorber lines
    - update_absorber_redshift(row, z_abs) - Update when redshift is changed manually
    - renumber_absorber_lines(old_row, new_row) - Handle row renumbering after deletion
    """
        
    def __init__(self, parent=None, colors=None):
        super().__init__(parent)
        self.parent = parent
        
        # Import color utility and use those colors if not provided
        self.colors = colors or list(clr.keys())[1:]  # Skip the first color (usually background)        

        # Define available line lists
        self.line_options = ['LLS', 'LLS Small', 'DLA', 'LBG', 'Gal', 'Eiger_Strong','AGN', 'None']
        
        # Initialize data storage for absorbers with IsPlotted state
        self.absorbers_df = pd.DataFrame(columns=['Zabs', 'LineList', 'Color', 'IsPlotted'])
        
        # Initialize UI components
        self.initUI()
        
        # Set up the dark theme after creating the widgets
        self.setup_dark_theme()
        
    def initUI(self):
        """Initialize the user interface components"""
        self.layout = QVBoxLayout(self)
        
        # Create table for displaying absorbers (reduced to 5 columns)
        self.table = QTableWidget(1, 5)  # Start with 1 row, will expand as needed
        self.table.setHorizontalHeaderLabels(['LineList', 'z', 'Color', 'Show', 'Remove'])
        
        # Configure table properties
        header = self.table.horizontalHeader()
        
        # Customize column widths to prevent overlap
        header.setSectionResizeMode(0, QHeaderView.Interactive)
        self.table.setColumnWidth(0, 80)  # Initial width
        header.setMinimumSectionSize(50)  # Prevents column from getting too small

        header.setSectionResizeMode(1, QHeaderView.Interactive)  # Redshift column
        header.setSectionResizeMode(2, QHeaderView.Interactive)  # Color column

        # Fixed width for action columns
        for col in range(3, 5):
            header.setSectionResizeMode(col, QHeaderView.Interactive)
            self.table.setColumnWidth(col, 30)  # Adjust width as needed
        
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
        line_combo.setStyleSheet("""
            QComboBox { 
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 6px;
                padding: 6px;
                font-size: 10px;
            }
            QComboBox::drop-down {
                background-color: #48484A;
                border-left: 1px solid #636366;
                width: 30px;
            }
            QComboBox QAbstractItemView {
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                selection-background-color: #0A84FF;
                selection-color: white;
            }
        """)
        self.table.setCellWidget(row_index, 0, line_combo)
        
        # Color combo box
        color_combo = QComboBox()
        for color in self.colors:
            color_combo.addItem(color)
        color_combo.setStyleSheet("""
            QComboBox { 
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 6px;
                padding: 6px;
                font-size: 14px;
            }
            QComboBox::drop-down {
                background-color: #48484A;
                border-left: 1px solid #636366;
                width: 30px;
            }
            QComboBox QAbstractItemView {
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                selection-background-color: #0A84FF;
                selection-color: white;
            }
        """)
        self.table.setCellWidget(row_index, 2, color_combo)
        
        # Show checkbox with direct row reference
        show_check = QCheckBox()
        show_check.setStyleSheet("""
            QCheckBox {
                color: #F2F2F7;
                spacing: 5px;
            }
            QCheckBox::indicator {
                width: 18px;
                height: 18px;
                border-radius: 4px;
            }
            QCheckBox::indicator:unchecked {
                background-color: #3A3A3C;
                border: 1px solid #636366;
            }
            QCheckBox::indicator:checked {
                background-color: #0A84FF;
                border: 1px solid #636366;
            }
        """)
        # Store the row index directly on the checkbox
        show_check.row_index = row_index
        show_check.stateChanged.connect(self.toggle_absorber_visibility_from_checkbox)
        self.table.setCellWidget(row_index, 3, show_check)
        
        # Remove button with direct row reference
        remove_btn = QPushButton("✕", self)
        remove_btn.setStyleSheet("""
            QPushButton {
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 6px;
                padding: 2px;
                font-size: 10px;
            }
            QPushButton:hover {
                background-color: #FF453A;
                color: white;
            }
        """)
        # Store the row index directly on the button
        remove_btn.row_index = row_index
        remove_btn.clicked.connect(self.remove_absorber_from_button)
        self.table.setCellWidget(row_index, 4, remove_btn)
    
    def toggle_absorber_visibility_from_checkbox(self, state):
        """Handle checkbox state change to plot or hide absorber"""
        checkbox = self.sender()
        if hasattr(checkbox, 'row_index'):
            row = checkbox.row_index
            is_checked = (state == QtCore.Qt.Checked)
            
            # Get the absorber data
            z_abs = float(self.table.item(row, 1).text()) if self.table.item(row, 1) else None
            if not z_abs:
                return
                
            line_list = self.table.cellWidget(row, 0).currentText()
            color = self.table.cellWidget(row, 2).currentText()
            
            # Check if this absorber has been plotted before
            is_plotted = False
            if row < len(self.absorbers_df):
                is_plotted = self.absorbers_df.at[row, 'IsPlotted'] if 'IsPlotted' in self.absorbers_df.columns else False
            
            if is_checked:
                # Show/plot the absorber
                if not is_plotted:
                    # If not previously plotted, plot it now
                    if hasattr(self.parent, 'plot_absorber_lines'):
                        success = self.parent.plot_absorber_lines(row, z_abs, line_list, color)
                        if success and row < len(self.absorbers_df):
                            self.absorbers_df.at[row, 'IsPlotted'] = True
                else:
                    # If already plotted but hidden, show it
                    if hasattr(self.parent, 'toggle_absorber_visibility'):
                        # Call the method with only the arguments it expects
                        self.parent.toggle_absorber_visibility(row)  # No third argument
            else:
                # Hide the absorber
                if is_plotted:
                    if hasattr(self.parent, 'toggle_absorber_visibility'):
                        # Call the method with only the arguments it expects
                        self.parent.toggle_absorber_visibility(row)  # No third argument

    
    def remove_absorber_from_button(self):
        """Remove absorber using row index stored on the button"""
        button = self.sender()
        if hasattr(button, 'row_index'):
            row = button.row_index
            print(f"Removing absorber at row {row}")
            self.remove_absorber(row)
        else:
            print("Button does not have a row_index property")
    
    def add_absorber(self, z_abs, line_list=None, color=None):
        """Add a new absorber to the manager"""
        # Default values if not provided
        line_list = line_list or 'LLS'
        color = color or 'white'
        
        # Add to dataframe with IsPlotted = True
        new_row = pd.Series({'Zabs': z_abs, 'LineList': line_list, 
                             'Color': color, 'IsPlotted': True})
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
            
            # Important: Add widgets BEFORE populating to ensure button connections exist
            self._add_row_widgets(new_row_idx)
            self._populate_row(new_row_idx, z_abs, line_list, color)
            
        
        # Always ensure there's at least one empty row at the end for future entries
        last_row = self.table.rowCount() - 1
        if self.table.item(last_row, 1) is not None:  # If the last row is not empty
            empty_row_idx = self.table.rowCount()
            self.table.setRowCount(empty_row_idx + 1)
            self._add_row_widgets(empty_row_idx)

            
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
    
    def remove_absorber(self, row):
        """Remove the absorber at the specified row"""
        try:
            # Call parent's method to remove absorber lines if available
            if hasattr(self.parent, 'remove_absorber_lines'):
                self.parent.remove_absorber_lines(row)
            
            # Need to get data before removing the row
            if row < len(self.absorbers_df):
                # Remove from dataframe
                self.absorbers_df = self.absorbers_df.drop(row).reset_index(drop=True)
            
            # Remove the row from the table
            self.table.removeRow(row)
            
            # Update row_index for all buttons after the removed row
            for r in range(row, self.table.rowCount()):
                for c in range(3, 5):  # Checkbox and remove button columns
                    widget = self.table.cellWidget(r, c)
                    if widget and hasattr(widget, 'row_index'):
                        # Update the row index to match the new position
                        widget.row_index = r
            
            # Add a new empty row at the end to maintain minimum rows
            new_row_idx = self.table.rowCount()
            self.table.insertRow(new_row_idx)
            self._add_row_widgets(new_row_idx)
                
            return True
        except Exception as e:
            print(f"Error removing absorber: {e}")
            return False
    
    def on_cell_changed(self, row, column):
        """Handle manual edits to cells"""
        if column == 1:  # Redshift column
            try:
                # Only proceed if there's actually text in the cell
                if self.table.item(row, 1) is not None and self.table.item(row, 1).text():
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
        Apply a dark theme matching the main application's color scheme
        """
        self.setStyleSheet("""
        QWidget {
            background-color: #353535;
            color: #F2F2F7;
        }
        
        QTableWidget {
            background-color: #252525;
            alternate-background-color: #353535;
            color: #F2F2F7;
            selection-background-color: #0A84FF;
            gridline-color: #636366;
            border-radius: 6px;
        }
        
        QTableWidget::item {
            background-color: #252525;
            color: #F2F2F7;
            padding: 4px;
        }
        
        QTableWidget::item:selected {
            background-color: #0A84FF;
            color: white;
        }
        
        QHeaderView::section {
            background-color: #3A3A3C;
            color: #F2F2F7;
            padding: 6px;
            border: 1px solid #636366;
            font-size: 14px;
        }
        
        QPushButton {
            background-color: #3A3A3C;
            color: #F2F2F7;
            border: 1px solid #636366;
            border-radius: 6px;
            padding: 6px;
            font-size: 14px;
        }
        
        QPushButton:hover {
            background-color: #48484A;
            color: white;
        }
        
        QCheckBox {
            color: #F2F2F7;
            spacing: 5px;
        }
        """)
        
        # Set alternating row colors for table
        self.table.setAlternatingRowColors(True)