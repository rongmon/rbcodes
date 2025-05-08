import pandas as pd
from PyQt5.QtWidgets import (QWidget, QTableWidget, QTableWidgetItem, 
                             QHBoxLayout, QVBoxLayout, QComboBox, QPushButton,
                             QHeaderView, QSizePolicy, QApplication, QCheckBox)
from PyQt5.QtGui import QColor
from PyQt5 import QtCore, QtGui, QtWidgets

class AbsorberManager(QWidget):
    """
    A widget for managing multiple absorber systems with their redshifts, line lists, and colors.
    Allows adding, removing, plotting, and hiding absorbers using a unified checkbox approach.
    
    Note: To fully support this widget, the parent should implement:
    - plot_absorber_lines(row, z_abs, line_list, color) - Plot lines for this absorber
    - remove_absorber_lines(row) - Remove lines for this absorber
    - toggle_absorber_visibility(row, visible) - Toggle or set visibility of absorber lines
    - update_absorber_redshift(row, z_abs) - Update when redshift is changed manually
    """
        
    def __init__(self, parent=None, colors=None):
        super().__init__(parent)
        self.parent = parent
        
        # Define theme colors
        self.theme = {
            'bg': '#252525',
            'widget_bg': '#353535',
            'header_bg': '#3A3A3C',
            'input_bg': '#3A3A3C',
            'border': '#636366',
            'text': '#F2F2F7',
            'accent': '#0A84FF',
            'danger': '#FF453A',
            'hover': '#48484A'
        }
        
        # Use provided colors or default list
        try:
            # Try to import color utility if available
            from rbcodes.utils import rb_utility as rt
            clr = rt.rb_set_color()
            self.colors = colors or list(clr.keys())[1:]  # Skip first color (background)
        except ImportError:
            # Fallback colors if utility not available
            self.colors = colors or ['white', 'red', 'blue', 'green', 'yellow', 'cyan', 'magenta']

        # Define available line lists
        self.line_options = ['LLS', 'LLS Small', 'DLA', 'LBG', 'Gal', 'Eiger_Strong', 'AGN', 'None']
        
        # Initialize data storage for absorbers
        self.absorbers_df = pd.DataFrame(columns=['Zabs', 'LineList', 'Color', 'IsPlotted'])
        
        # Initialize UI components
        self.initUI()
    
    def initUI(self):
        """Initialize the user interface with optimized layout"""
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(5, 5, 5, 5)
        
        # Create table for displaying absorbers
        self.table = QTableWidget(1, 5)
        self.table.setHorizontalHeaderLabels(['LineList', 'z', 'Color', 'Show', 'Remove'])
        
        # Set column properties
        header = self.table.horizontalHeader()
        column_widths = [80, 70, 70, 30, 30]
        
        for i, width in enumerate(column_widths):
            header.setSectionResizeMode(i, QHeaderView.Interactive)
            self.table.setColumnWidth(i, width)
        
        header.setMinimumSectionSize(30)
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.table.setAlternatingRowColors(True)
        
        # Add initial rows with widgets
        self._add_row(0)
        self._add_row(1)
        
        # Add table to layout
        self.layout.addWidget(self.table)
        
        # Connect signals
        self.table.cellChanged.connect(self.on_cell_changed)
        
        # Apply theme
        self.apply_theme()
    
    def _add_row(self, row_index):
        """Add a new row with all necessary widgets"""
        # Create widgets
        line_combo = self._create_combo_box(self.line_options)
        color_combo = self._create_combo_box(self.colors)
        show_check = self._create_checkbox()
        remove_btn = self._create_button("✕")
        
        # Connect signals with lambdas to capture row context
        show_check.stateChanged.connect(
            lambda state, row=row_index: self.toggle_absorber_visibility_from_checkbox(row, state)
        )
        remove_btn.clicked.connect(
            lambda _, row=row_index: self.remove_absorber(row)
        )
        
        # Add widgets to table
        self.table.setCellWidget(row_index, 0, line_combo)
        self.table.setCellWidget(row_index, 2, color_combo)
        self.table.setCellWidget(row_index, 3, show_check)
        self.table.setCellWidget(row_index, 4, remove_btn)
    
    def _create_combo_box(self, items):
        """Factory method for creating styled combo boxes"""
        combo = QComboBox()
        for item in items:
            combo.addItem(item)
        return combo
    
    def _create_checkbox(self):
        """Factory method for creating styled checkboxes"""
        return QCheckBox()
    
    def _create_button(self, text):
        """Factory method for creating styled buttons"""
        btn = QPushButton(text, self)
        # Make remove button's text color stand out for clarity
        if text == "✕":
            btn.setProperty("type", "danger")
        return btn
    
    def toggle_absorber_visibility_from_checkbox(self, row, state):
        """Handle checkbox state change to plot or hide absorber"""
        # Determine if checked
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
                    self.parent.toggle_absorber_visibility(row)
        else:
            # Hide the absorber
            if is_plotted and hasattr(self.parent, 'toggle_absorber_visibility'):
                self.parent.toggle_absorber_visibility(row)
    
    def add_absorber(self, z_abs, line_list=None, color=None):
        """Add a new absorber to the manager"""
        # Default values if not provided
        line_list = line_list or 'LLS'
        color = color or 'white'
        
        # Add to dataframe with IsPlotted = True
        new_row = pd.Series({'Zabs': z_abs, 'LineList': line_list, 
                            'Color': color, 'IsPlotted': True})
        self.absorbers_df = pd.concat([self.absorbers_df, new_row.to_frame().T], 
                                      ignore_index=True)
        
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
            
            # Add widgets before populating to ensure connections exist
            self._add_row(new_row_idx)
            self._populate_row(new_row_idx, z_abs, line_list, color)
        
        # Ensure there's at least one empty row at the end
        last_row = self.table.rowCount() - 1
        if self.table.item(last_row, 1) is not None:  # Last row not empty
            empty_row_idx = self.table.rowCount()
            self.table.setRowCount(empty_row_idx + 1)
            self._add_row(empty_row_idx)
            
        # Return success
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
            
            # Remove from dataframe if row exists
            if row < len(self.absorbers_df):
                self.absorbers_df = self.absorbers_df.drop(row).reset_index(drop=True)
            
            # Remove the row from the table
            self.table.removeRow(row)
            
            # Add a new empty row at the end if needed
            if self.table.rowCount() == 0 or (
                self.table.rowCount() > 0 and
                self.table.item(self.table.rowCount()-1, 1) is not None
            ):
                new_row_idx = self.table.rowCount()
                self.table.insertRow(new_row_idx)
                self._add_row(new_row_idx)
                
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
    
    def apply_theme(self):
        """Apply centralized theme to all components"""
        t = self.theme  # Short alias for theme colors
        
        # Base widget style
        widget_style = f"""
            QWidget {{
                background-color: {t['widget_bg']};
                color: {t['text']};
            }}
        """
        
        # Table style
        table_style = f"""
            QTableWidget {{
                background-color: {t['bg']};
                alternate-background-color: {t['widget_bg']};
                color: {t['text']};
                selection-background-color: {t['accent']};
                gridline-color: {t['border']};
                border-radius: 6px;
            }}
            
            QTableWidget::item {{
                background-color: {t['bg']};
                color: {t['text']};
                padding: 4px;
            }}
            
            QTableWidget::item:selected {{
                background-color: {t['accent']};
                color: white;
            }}
            
            QHeaderView::section {{
                background-color: {t['header_bg']};
                color: {t['text']};
                padding: 6px;
                border: 1px solid {t['border']};
                font-size: 10px;
            }}
        """
        
        # Combobox style
        combobox_style = f"""
            QComboBox {{ 
                background-color: {t['input_bg']};
                color: {t['text']};
                border: 1px solid {t['border']};
                border-radius: 6px;
                padding: 6px;
                font-size: 10px;
            }}
            QComboBox::drop-down {{
                background-color: {t['hover']};
                border-left: 1px solid {t['border']};
                width: 30px;
            }}
            QComboBox QAbstractItemView {{
                background-color: {t['input_bg']};
                color: {t['text']};
                border: 1px solid {t['border']};
                selection-background-color: {t['accent']};
                selection-color: white;
            }}
        """
        
        # Button style
        button_style = f"""
            QPushButton {{
                background-color: {t['input_bg']};
                color: {t['text']};
                border: 1px solid {t['border']};
                border-radius: 6px;
                padding: 2px;
                font-size: 10px;
            }}
            QPushButton:hover {{
                background-color: {t['hover']};
                color: white;
            }}
            QPushButton[type="danger"]:hover {{
                background-color: {t['danger']};
                color: white;
            }}
        """
        
        # Checkbox style
        checkbox_style = f"""
            QCheckBox {{
                color: {t['text']};
                spacing: 5px;
            }}
            QCheckBox::indicator {{
                width: 18px;
                height: 18px;
                border-radius: 4px;
            }}
            QCheckBox::indicator:unchecked {{
                background-color: {t['input_bg']};
                border: 1px solid {t['border']};
            }}
            QCheckBox::indicator:checked {{
                background-color: {t['accent']};
                border: 1px solid {t['border']};
            }}
        """
        
        # Combine all styles
        combined_style = "\n".join([
            widget_style, 
            table_style, 
            combobox_style, 
            button_style, 
            checkbox_style
        ])
        
        # Apply combined style to the widget
        self.setStyleSheet(combined_style)
        
        # Apply styles to comboboxes
        for row in range(self.table.rowCount()):
            for col in [0, 2]:  # LineList and Color combo columns
                widget = self.table.cellWidget(row, col)
                if isinstance(widget, QComboBox):
                    widget.setStyleSheet(combobox_style)