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
    Allows adding, removing, plotting, and hiding absorbers with a single checkbox.
    
    Note: To fully support this widget, the parent should implement:
    - plot_absorber_lines(row, z_abs, line_list, color) - Plot lines for this absorber
    - remove_absorber_lines(row) - Remove lines for this absorber
    - toggle_absorber_visibility(row, is_visible) - Toggle visibility of absorber lines
    - update_absorber_redshift(row, z_abs) - Update when redshift is changed manually
    """
        
    def __init__(self, parent=None, colors=None, absorbers_df=None):
        super().__init__(parent)
        self.parent = parent
        
        # Import color utility and use those colors if not provided
        self.colors = colors or list(clr.keys())[1:]  # Skip the first color (usually background)        

        # Define available line lists
        self.line_options = ['LLS', 'LLS Small', 'DLA', 'LBG', 'Gal', 'Eiger_Strong','AGN', 'None']
        
        # Initialize data storage for absorbers
        if absorbers_df is not None:
            # If DataFrame is provided, use it (must have Zabs, LineList, Color columns)
            self.absorbers_df = absorbers_df.copy()
            # Add Visible column if it doesn't exist
            if 'Visible' not in self.absorbers_df.columns:
                self.absorbers_df['Visible'] = False
        else:
            # Initialize empty DataFrame
            self.absorbers_df = pd.DataFrame(columns=['Zabs', 'LineList', 'Color', 'Visible'])
        
        # Initialize UI components
        self.initUI()
        
        # Set up the dark theme after creating the widgets
        self.setup_dark_theme()
        
        # If absorbers were provided, populate the table
        if absorbers_df is not None and not absorbers_df.empty:
            self.populate_table_from_df()
    
    def populate_table_from_df(self):
        """Populate the table from an existing DataFrame"""
        # Clear existing rows
        self.table.setRowCount(0)
        
        # Add rows for each absorber in the DataFrame
        for idx, row in self.absorbers_df.iterrows():
            row_idx = self.table.rowCount()
            self.table.insertRow(row_idx)
            self._add_row_widgets(row_idx)
            
            # Set values from DataFrame
            self._populate_row(
                row_idx, 
                row['Zabs'], 
                row['LineList'], 
                row['Color']
            )
            
            # Set checkbox state if Visible column exists
            if 'Visible' in self.absorbers_df.columns:
                visible = row['Visible']
                checkbox = self.table.cellWidget(row_idx, 3)
                if checkbox:
                    checkbox.setChecked(visible)
                    # Plot if visible is True
                    if visible:
                        self.plot_absorber(row_idx)
        
        # Add one empty row at the end
        empty_row_idx = self.table.rowCount()
        self.table.insertRow(empty_row_idx)
        self._add_row_widgets(empty_row_idx)
        
    def initUI(self):
        """Initialize the user interface components"""
        self.layout = QVBoxLayout(self)
        
        # Create table for displaying absorbers
        self.table = QTableWidget(1, 5)  # Start with 1 row, 5 columns
        self.table.setHorizontalHeaderLabels(['LineList', 'z', 'Color', 'Plot/Hide', 'Remove'])
        
        # Configure table properties
        header = self.table.horizontalHeader()
        
        # Customize column widths
        header.setSectionResizeMode(0, QHeaderView.Interactive)
        self.table.setColumnWidth(0, 80)  # Line Lists column
        header.setMinimumSectionSize(50)
        
        header.setSectionResizeMode(1, QHeaderView.Interactive)  # Redshift column
        header.setSectionResizeMode(2, QHeaderView.Interactive)  # Color column
        header.setSectionResizeMode(3, QHeaderView.Interactive)  # Plot/Hide checkbox column
        self.table.setColumnWidth(3, 60)
        
        # Fixed width for action buttons
        header.setSectionResizeMode(4, QHeaderView.Interactive)
        self.table.setColumnWidth(4, 60)  # Remove button
        
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Set default widgets for initial row
        self._add_row_widgets(0)
        
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
                border: 0.5px solid #2A82DA;
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
                border: 0.5px solid #2A82DA;
            }
            QComboBox::drop-down {
                background-color: #2A82DA;
            }
        """)
        self.table.setCellWidget(row_index, 2, color_combo)
        
        # Plot/Hide checkbox
        checkbox = QCheckBox()
        checkbox.setStyleSheet("""
            QCheckBox {
                background-color: transparent;
                color: white;
            }
            QCheckBox::indicator {
                width: 12px;
                height: 12px;
            }
            QCheckBox::indicator:unchecked {
                border: 0.5px solid #2A82DA;
                background-color: #353535;
            }
            QCheckBox::indicator:checked {
                background-color: #2A82DA;
                border: 1px solid #2A82DA;
            }
        """)
        # Store the row index directly on the checkbox
        checkbox.row_index = row_index
        checkbox.stateChanged.connect(self.toggle_absorber_from_checkbox)
        
        # Create a container widget for centering the checkbox
        checkbox_container = QWidget()
        checkbox_layout = QHBoxLayout(checkbox_container)
        checkbox_layout.addWidget(checkbox)
        checkbox_layout.setAlignment(QtCore.Qt.AlignCenter)
        checkbox_layout.setContentsMargins(0, 0, 0, 0)
        
        self.table.setCellWidget(row_index, 3, checkbox_container)
        
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
        self.table.setCellWidget(row_index, 4, remove_btn)
    
    def toggle_absorber_from_checkbox(self, state):
        """Handle checkbox state change"""
        checkbox = self.sender()
        if hasattr(checkbox, 'row_index'):
            row = checkbox.row_index
            if self.table.item(row, 1) is None:  # Skip if row is empty
                return
                
            is_checked = (state == QtCore.Qt.Checked)
            
            # Update the absorbers_df
            if row < len(self.absorbers_df):
                self.absorbers_df.at[row, 'Visible'] = is_checked
            
            if is_checked:
                # Plot the absorber
                self.plot_absorber(row)
            else:
                # Call parent's method to hide absorber
                if hasattr(self.parent, 'remove_absorber_lines'):
                    self.parent.remove_absorber_lines(row)
    
    def remove_absorber_from_button(self):
        """Remove absorber using row index stored on the button"""
        button = self.sender()
        if hasattr(button, 'row_index'):
            row = button.row_index
            print(f"Removing absorber at row {row}")
            self.remove_absorber(row)
        else:
            print("Button does not have a row_index property")
    
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
    
    def add_absorber(self, z_abs, line_list=None, color=None, visible=False):
        """Add a new absorber to the manager"""
        # Default values if not provided
        line_list = line_list or 'LLS'
        color = color or 'white'
        
        # Add to dataframe
        new_row = pd.Series({
            'Zabs': z_abs, 
            'LineList': line_list, 
            'Color': color,
            'Visible': visible
        })
        self.absorbers_df = self.absorbers_df.append(new_row, ignore_index=True)
        
        # Find the first empty row or add a new one
        row_found = False
        for row in range(self.table.rowCount()):
            if self.table.item(row, 1) is None:
                # This is an empty row, use it
                self._populate_row(row, z_abs, line_list, color)
                # Set checkbox state
                checkbox_container = self.table.cellWidget(row, 3)
                if checkbox_container:
                    for child in checkbox_container.children():
                        if isinstance(child, QCheckBox):
                            child.setChecked(visible)
                            break
                row_found = True
                break
        
        if not row_found:
            # No empty rows, add a new one
            new_row_idx = self.table.rowCount()
            self.table.setRowCount(new_row_idx + 1)
            
            # Important: Add widgets BEFORE populating
            self._add_row_widgets(new_row_idx)
            self._populate_row(new_row_idx, z_abs, line_list, color)
            
            # Set checkbox state
            checkbox_container = self.table.cellWidget(new_row_idx, 3)
            if checkbox_container:
                for child in checkbox_container.children():
                    if isinstance(child, QCheckBox):
                        child.setChecked(visible)
                        break
            
        # Always ensure there's at least one empty row at the end
        last_row = self.table.rowCount() - 1
        if self.table.item(last_row, 1) is not None:  # If the last row is not empty
            empty_row_idx = self.table.rowCount()
            self.table.setRowCount(empty_row_idx + 1)
            self._add_row_widgets(empty_row_idx)
        
        # If visible is True, plot the absorber
        if visible:
            for row in range(self.table.rowCount()):
                if self.table.item(row, 1) is not None and float(self.table.item(row, 1).text()) == z_abs:
                    self.plot_absorber(row)
                    break
        
        # Return success status
        return True
    
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
                success = self.parent.plot_absorber_lines(row, z_abs, line_list, color, alpha=0.35, linewidth=0.5)
                
                # Update visibility in dataframe if successful
                if success and row < len(self.absorbers_df):
                    self.absorbers_df.at[row, 'Visible'] = True
                    
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
            
            # Update row_index for all buttons and checkboxes after the removed row
            for r in range(row, self.table.rowCount()):
                # Update remove button row_index
                btn = self.table.cellWidget(r, 4)
                if btn and hasattr(btn, 'row_index'):
                    btn.row_index = r
                
                # Update checkbox row_index
                checkbox_container = self.table.cellWidget(r, 3)
                if checkbox_container:
                    for child in checkbox_container.children():
                        if isinstance(child, QCheckBox) and hasattr(child, 'row_index'):
                            child.row_index = r
                            break
            
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
                    else:
                        # New row data
                        line_list = self.table.cellWidget(row, 0).currentText()
                        color = self.table.cellWidget(row, 2).currentText()
                        
                        # Check if checkbox is checked
                        visible = False
                        checkbox_container = self.table.cellWidget(row, 3)
                        if checkbox_container:
                            for child in checkbox_container.children():
                                if isinstance(child, QCheckBox):
                                    visible = child.isChecked()
                                    break
                        
                        # Add to dataframe
                        new_row = pd.Series({
                            'Zabs': z_abs, 
                            'LineList': line_list, 
                            'Color': color,
                            'Visible': visible
                        })
                        self.absorbers_df = self.absorbers_df.append(new_row, ignore_index=True)
                    
                    # Notify parent if needed
                    if hasattr(self.parent, 'update_absorber_redshift'):
                        self.parent.update_absorber_redshift(row, z_abs)
                    
                    # If checkbox is checked, replot the absorber
                    checkbox_container = self.table.cellWidget(row, 3)
                    if checkbox_container:
                        for child in checkbox_container.children():
                            if isinstance(child, QCheckBox) and child.isChecked():
                                self.plot_absorber(row)
                                break
                    
                    # Ensure there's an empty row at the end
                    last_row = self.table.rowCount() - 1
                    if self.table.item(last_row, 1) is not None:
                        empty_row_idx = self.table.rowCount()
                        self.table.setRowCount(empty_row_idx + 1)
                        self._add_row_widgets(empty_row_idx)
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
    
    def get_all_absorber_data(self):
        """Get all absorber data as a DataFrame"""
        return self.absorbers_df.copy()

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
        
        QCheckBox::indicator:checked {
            background-color: #2A82DA;
        }
        
        QCheckBox::indicator:unchecked {
            background-color: #353535;
            border: 1px solid #636366;
        }
        """)