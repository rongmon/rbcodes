import sys
import pandas as pd

from PyQt5.QtWidgets import (QWidget, QTextEdit, QVBoxLayout, QApplication, QSizePolicy, 
                             QLabel, QHBoxLayout)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor, QPalette, QFont

class MessageBox(QWidget):
    """
    A widget that displays messages to the user with support for colored text.
    Designed to match the style of PlotSpec_Integrated.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Dark theme colors
        self.dark_background = '#353535'      # Main dark background
        self.darker_background = '#252525'    # Slightly darker background for inputs
        self.white_text = '#FFFFFF'           # White text
        self.highlight_color = '#2A82DA'      # Highlight/accent color
        
        # Configure the widget
        self.initUI()
        
    def initUI(self):
        # Set size policy
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setMinimumWidth(300)
        self.setMinimumHeight(100)
        
        # Create layout
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Add a title label
        title_layout = QHBoxLayout()
        title_label = QLabel("Message Window")
        title_label.setStyleSheet(f"""
            background-color: {self.darker_background}; 
            color: {self.white_text};
            font-weight: bold;
            padding: 4px;
        """)
        title_layout.addWidget(title_label)
        
        # Create the text edit widget for messages
        self.text_edit = QTextEdit()
        self.text_edit.setReadOnly(True)
        self.text_edit.setPlaceholderText("This message box will display important messages.")
        
        # Set the default text with a welcome message
        self.text_edit.setHtml('<span style="color:#FFFFFF;">This message box will display important messages.</span>')
        
        # Apply stylesheet to match PlotSpec_Integrated style
        self.setStyleSheet(f"""
            QWidget {{
                background-color: {self.dark_background};
            }}
            
            QTextEdit {{
                background-color: {self.darker_background};
                color: {self.white_text};
                border: 1px solid {self.highlight_color};
                border-radius: 3px;
                padding: 5px;
                font-size: 12px;
                selection-background-color: {self.highlight_color};
            }}
        """)
        
        # Add widgets to layout
        layout.addLayout(title_layout)
        layout.addWidget(self.text_edit)
        
        # Default message
        self.message = ''
        
    def on_sent_message(self, sent_message, hexColor='#FFFFFF'):
        """
        Send a message with optional color
        
        :param sent_message: Message text
        :param hexColor: Hex color code for the message (default white)
        """
        # Sanitize the message to prevent HTML injection
        sanitized_message = sent_message.replace('<', '&lt;').replace('>', '&gt;')
        
        # Create a styled HTML message
        message = f'<span style="color:{hexColor};">{sanitized_message}</span>'
        
        # Set the message with HTML formatting
        self.text_edit.setHtml(message)
        
        # Optional: scroll to the bottom to ensure the message is visible
        self.text_edit.verticalScrollBar().setValue(
            self.text_edit.verticalScrollBar().maximum()
        )
        
    def append_message(self, sent_message, hexColor='#FFFFFF'):
        """
        Append a message with optional color to the existing messages
        
        :param sent_message: Message text
        :param hexColor: Hex color code for the message (default white)
        """
        # Sanitize the message to prevent HTML injection
        sanitized_message = sent_message.replace('<', '&lt;').replace('>', '&gt;')
        
        # Create a styled HTML message
        message = f'<span style="color:{hexColor};">{sanitized_message}</span><br>'
        
        # Append the message with HTML formatting
        current_html = self.text_edit.toHtml()
        
        # Insert the new message before the closing body and html tags
        if '</body>' in current_html:
            new_html = current_html.replace('</body>', message + '</body>')
            self.text_edit.setHtml(new_html)
        else:
            self.text_edit.append(message)
        
        # Scroll to the bottom
        self.text_edit.verticalScrollBar().setValue(
            self.text_edit.verticalScrollBar().maximum()
        )
        
    def clear(self):
        """Clear all messages from the text edit widget"""
        self.text_edit.clear()


# For testing the widget in isolation
if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # Apply dark theme to entire application
    app.setStyle("Fusion")
    dark_palette = QPalette()
    dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.WindowText, Qt.white)
    dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
    dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
    dark_palette.setColor(QPalette.ToolTipText, Qt.white)
    dark_palette.setColor(QPalette.Text, Qt.white)
    dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ButtonText, Qt.white)
    dark_palette.setColor(QPalette.BrightText, Qt.red)
    dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(dark_palette)
    
    # Create widget and show it
    message_box = MessageBox()
    message_box.show()
    
    # Show some test messages
    message_box.on_sent_message("This is a test message", "#FFFFFF")
    message_box.append_message("This is a success message", "#00FF00")
    message_box.append_message("This is an error message", "#FF0000")
    message_box.append_message("This is a warning message", "#FFA500")
    
    sys.exit(app.exec_())