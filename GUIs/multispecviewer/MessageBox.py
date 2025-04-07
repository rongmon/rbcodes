import sys
import pandas as pd

from PyQt5.QtWidgets import QWidget, QTextEdit, QVBoxLayout
from PyQt5.QtCore import pyqtSignal

from PyQt5 import QtCore
from PyQt5 import QtGui

class MessageBox(QWidget):
    def __init__(self):
        super().__init__()

        # Dark theme colors
        dark_background = '#353535'      # Main dark background
        darker_background = '#252525'    # Slightly darker background for inputs
        white_text = '#FFFFFF'           # White text
        highlight_color = '#2A82DA'      # Highlight/accent color

        # Dark theme stylesheet
        # Fix around line 14-33
        dark_stylesheet = """
        QWidget {
            background-color: """ + dark_background + """;
            color: """ + white_text + """;
        }
        
        QTextEdit {
            background-color: """ + darker_background + """;
            color: """ + white_text + """;
            border: 1px solid """ + highlight_color + """;
            padding: 5px;
        }
        """

        self.message = ''

        self.te = QTextEdit()
        self.te.setPlaceholderText('This message box will display important messages.')
        self.te.setReadOnly(True)

        # Rich text support for colored messages
        self.te.setHtml('<span style="color:#FFFFFF;">This message box will display important messages.</span>')

        self.setFixedHeight(200)
        self.setFixedWidth(250)

        layout = QVBoxLayout()
        layout.addWidget(self.te)
        layout.setAlignment(QtCore.Qt.AlignCenter)

        # Apply dark theme stylesheet
        self.setStyleSheet(dark_stylesheet)

        self.setLayout(layout)

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
        self.te.setHtml(message)
        
        # Optional: scroll to the bottom
        self.te.verticalScrollBar().setValue(
            self.te.verticalScrollBar().maximum()
        )