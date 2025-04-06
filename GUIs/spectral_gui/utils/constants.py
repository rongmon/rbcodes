#!/usr/bin/env python
"""
Constants - Shared constants for the spectrum viewer application.
"""

# Speed of light in km/s
SPEED_OF_LIGHT = 2.9979e5

# Available line lists
LINE_LISTS = [
    'LLS',
    'LLS Small',
    'DLA',
    'LBG',
    'Gal',
    'None'
]

# Color definitions - maps color names to RGB values
# Create a colors dictionary matching the original
COLORS = {
    'white': '#FFFFFF',
    'yellow': '#FFFF00',
    'orange': '#FFA500',
    'red': '#FF0000',
    'pale_red': '#FF9999',
    'green': '#00FF00',
    'blue': '#0000FF',
    'cyan': '#00FFFF',
    'magenta': '#FF00FF',
    'pink': '#FFC0CB',
    'purple': '#800080',
    'teal': '#008080',
    'light_gray': '#D3D3D3',
    'gray': '#808080',
    'black': '#000000'
}

# Default wavelength range for plotting
DEFAULT_WAVE_RANGE = [3000, 9000]  # Angstroms

# Default velocity range for line inspection
DEFAULT_VELOCITY_RANGE = [-1000, 1000]  # km/s

# Default plot settings
DEFAULT_PLOT_SETTINGS = {
    'linewidth': 0.9,
    'marker_size': 5,
    'font_size': 10,
    'error_alpha': 0.7
}

# Status flags for line identification
LINE_STATUS = {
    0: 'Non-Detection',
    1: 'Detection',
    2: 'Blended Detection',
    3: 'Low-Confidence Detection'
}

# Line status display settings
LINE_STATUS_STYLES = {
    0: {'color': 'light_gray', 'linestyle': '--'},
    1: {'color': 'yellow', 'linestyle': '--'},
    2: {'color': 'cyan', 'linestyle': '-.'},
    3: {'color': 'pale_red', 'linestyle': ':'}
}

# Special ion groups
ION_GROUPS = {
    'HI': ['HI 1215', 'HI 1025', 'HI 972', 'HI 949', 'HI 937'],
    'CIV': ['CIV 1548', 'CIV 1550'],
    'SiIV': ['SiIV 1393', 'SiIV 1402'],
    'OVI': ['OVI 1031', 'OVI 1037'],
    'MgII': ['MgII 2796', 'MgII 2803'],
    'FeII': ['FeII 2600', 'FeII 2586', 'FeII 2382', 'FeII 2374', 'FeII 2344']
}

# Help text sections - for context-sensitive help
HELP_SECTIONS = {
    'main': """
        Main Window Help:
        Use keyboard shortcuts for most operations. See full help with 'H'.
    """,
    
    'absorber_manager': """
        Absorber Manager Help:
        Add, edit, and manage absorption systems.
        Enter a redshift in the last row to create a new absorber.
    """,
    
    'vstack': """
        VStack Display Help:
        Navigate with '<' and '>' keys.
        Toggle detection status with 'w' key.
        Press 'S' to save and exit.
    """,
    
    'measurements': """
        Measurement Help:
        Equivalent Width: Press 'E' at two points to measure.
        Gaussian Fit: Press 'G' at three points to fit.
    """
}
