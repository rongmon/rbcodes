"""
Shared pytest configuration for rbcodes.

GUI modules are imported and exercised without a display: Qt uses the
offscreen platform and matplotlib the Agg backend. Both must be set before
PyQt5 or pyplot is first imported.
"""
import os

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
os.environ.setdefault('MPLBACKEND', 'Agg')
