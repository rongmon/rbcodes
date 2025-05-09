import pandas as pd
from PyQt5.QtWidgets import (QWidget, QTableWidget, QTableWidgetItem, 
                             QHBoxLayout, QVBoxLayout, QComboBox, QPushButton,
                             QHeaderView, QSizePolicy, QApplication)
from PyQt5.QtGui import QColor
from PyQt5 import QtCore, QtGui, QtWidgets
from rbcodes.utils import rb_utility as rt
clr = rt.rb_set_color()

class AbsorberManager(QWidget):
    """
    A widget for managing multiple absorber systems with their redshifts, line lists, and colors.
    Allows adding, removing, plotting, and hiding absorbers.
    
    Note: To fully support this widget, the parent should implement:
    - plot_absorber_lines(row, z_abs, line_list, color) - Plot lines for this absorber
    - remove_absorber_lines(row) - Remove lines for this absorber
    - toggle_absorber_visibility(row) - Toggle visibility of absorber lines
    - update_absorber_redshift(row, z_abs) - Update when redshift is changed manually
    - renumber_absorber_lines(old_row, new_row) - Handle row renumbering after deletion
    """
        
    def __init__(self, parent=None, colors=None):
        super().__init__(parent)
        # Set up the dark theme after creating the widgets
        self.setup_dark_theme()
        self.parent = parent
        
        # Import color utility and use those colors if not provided
        self.colors = colors or list(clr.keys())[1:]  # Skip the first color (usually background)        

        # Define available line lists
        self.line_options = ['LLS', 'LLS Small', 'DLA', 'LBG', 'Gal', 'Eiger_Strong','AGN', 'None']
        # Initialize data storage for absorbers
        self.absorbers_df = pd.DataFrame(columns=['Zabs', 'LineList', 'Color'])
        
        # Initialize UI components
        self.initUI()
        
    def initUI(self):
        """Initialize the user interface components"""
        self.layout = QVBoxLayout(self)
        
        # Create table for displaying absorbers
        self.table = QTableWidget(1, 6)  # Start with 1 row, will expand as needed
        self.table.setHorizontalHeaderLabels(['LineList', 'z', 'Color', 'Plot', 'Hide','Remove'])
        
        # Configure table properties
        header = self.table.horizontalHeader()
        
        # Customize column widths to prevent overlap
        #header.setSectionResizeMode(0, QHeaderView.Stretch)  # Line Lists column stretches
        header.setSectionResizeMode(0, QHeaderView.Interactive)
        self.table.setColumnWidth(0, 80)  # Initial width
        header.setMinimumSectionSize(50)  # Prevents column from getting too small

        header.setSectionResizeMode(1, QHeaderView.Interactive)  # Redshift column
        header.setSectionResizeMode(2, QHeaderView.Interactive)  # Color column

        
        # Fixed width for action buttons
        for col in range(3, 6):
            header.setSectionResizeMode(col, QHeaderView.Interactive)
            self.table.setColumnWidth(col, 30)  # Adjust width as needed
        
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Set default widgets for initial rows
        self._add_row_widgets(0)
        #self._add_row_widgets(1)
        
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
                background-color: #353535; 
                color: white; 
                border: 1px solid #2A82DA;
            }
            QComboBox::drop-down {
                background-color: #2A82DA;
            }
        """)
        self.table.setCellWidget(row_index, 0, line_combo)
        
        # Color combo box
        color_combo = QComboBox()
        for color in self.colors:
            color_combo.addItem(color)
        color_combo.setStyleSheet("""
            QComboBox { 
                background-color: #353535; 
                color: white; 
                border: 1px solid #2A82DA;
            }
            QComboBox::drop-down {
                background-color: #2A82DA;
            }
        """)
        self.table.setCellWidget(row_index, 2, color_combo)
        
        # Plot button with direct row reference
        plot_btn = QPushButton("Plot", self)
        plot_btn.setStyleSheet("""
            QPushButton {
                background-color: #2A82DA; 
                color: white; 
                border: none;
                padding: 5px;
            }
            QPushButton:hover {
                background-color: #4A9FEA;
            }
        """)
        # Store the row index directly on the button
        plot_btn.row_index = row_index
        plot_btn.clicked.connect(self.plot_absorber_from_button)
        self.table.setCellWidget(row_index, 3, plot_btn)
        
        # Remove button with direct row reference
        remove_btn = QPushButton("Remove", self)
        remove_btn.setStyleSheet("""
            QPushButton {
                background-color: #FF6B6B; 
                color: white; 
                border: none;
                padding: 5px;
            }
            QPushButton:hover {
                background-color: #FF8A8A;
            }
        """)
        # Store the row index directly on the button
        remove_btn.row_index = row_index
        remove_btn.clicked.connect(self.remove_absorber_from_button)
        self.table.setCellWidget(row_index, 5, remove_btn)
        
        # Hide button with direct row reference
        hide_btn = QPushButton("Hide", self)
        hide_btn.setStyleSheet("""
            QPushButton {
                background-color: #FFA500; 
                color: white; 
                border: none;
                padding: 5px;
            }
            QPushButton:hover {
                background-color: #FFB740;
            }
        """)
        # Store the row index directly on the button
        hide_btn.row_index = row_index
        hide_btn.clicked.connect(self.toggle_hide_absorber_from_button)
        self.table.setCellWidget(row_index, 4, hide_btn)
    
    def get_index(self):
        """
        Gets the row and column of the clicked button.
        This is a comprehensive approach based on PlotSpec_Integrated.
        """
        # Get the button that was clicked
        button = QtWidgets.QApplication.focusWidget()
        if not button:
            return None, None
            
        # Try direct approach - check if the widget is in our table
        for row in range(self.table.rowCount()):
            for col in range(3, 6):  # Columns with buttons
                if self.table.cellWidget(row, col) == button:
                    print(f"Direct check - Button at row {row}, col {col}")
                    return row, col
        
        # Try indexAt with button position (might not work reliably)
        index = self.table.indexAt(button.pos())
        if index.isValid():
            row = index.row()
            col = index.column()
            print(f"IndexAt check - Button at row {row}, col {col}")
            return row, col
        
        print("Could not determine button position")
        return None, None
    
    def plot_absorber_from_button(self):
        """Plot absorber using row index stored on the button"""
        button = self.sender()
        if hasattr(button, 'row_index'):
            row = button.row_index
            print(f"Plotting absorber at row {row}")
            self.plot_absorber(row)
        else:
            print("Button does not have a row_index property")
    
    def remove_absorber_from_button(self):
        """Remove absorber using row index stored on the button"""
        button = self.sender()
        if hasattr(button, 'row_index'):
            row = button.row_index
            print(f"Removing absorber at row {row}")
            self.remove_absorber(row)
        else:
            print("Button does not have a row_index property")
    
    def toggle_hide_absorber_from_button(self):
        """Toggle hide absorber using row index stored on the button"""
        button = self.sender()
        if hasattr(button, 'row_index'):
            row = button.row_index
            print(f"Toggling visibility of absorber at row {row}")
            self.toggle_hide_absorber(row)
        else:
            print("Button does not have a row_index property")
        
    def get_button_row_index(self):
        """
        Gets the row index of the button that triggered the event.
        This technique is adapted from PlotSpec_Integrated.
        """
        button = self.sender()
        if not button:
            return -1
            
        # First try to find the button directly in the table
        for row in range(self.table.rowCount()):
            for col in range(3, 6):  # Columns 3, 4, 5 contain our buttons
                if self.table.cellWidget(row, col) == button:
                    return row
                    
        # Fallback: try using position
        index = self.table.indexAt(button.pos())
        if index.isValid():
            return index.row()
            
        return -1
    
    def plot_clicked(self):
        """Handle Plot button click"""
        row = self.get_button_row_index()
        if row >= 0:
            self.plot_absorber(row)
            
    def remove_clicked(self):
        """Handle Remove button click"""
        row = self.get_button_row_index()
        if row >= 0:
            self.remove_absorber(row)
            
    def hide_clicked(self):
        """Handle Hide button click"""
        row = self.get_button_row_index()
        if row >= 0:
            self.toggle_hide_absorber(row)
    
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
                success = self.parent.plot_absorber_lines(row, z_abs, line_list, color,alpha=0.35,linewidth=0.5)
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
                for c in range(3, 6):  # Button columns
                    btn = self.table.cellWidget(r, c)
                    if btn and hasattr(btn, 'row_index'):
                        # Update the row index to match the new position
                        btn.row_index = r
            
            # Add a new empty row at the end to maintain minimum rows
            new_row_idx = self.table.rowCount()
            self.table.insertRow(new_row_idx)
            self._add_row_widgets(new_row_idx)
                
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
        
        