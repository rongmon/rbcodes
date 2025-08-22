"""
Text utilities module for the absorption line analysis toolbox.
This module handles text formatting and display for measurements.
"""

import numpy as np

def plot_text(parent, line):
    """
    Generate formatted text for displaying measurement results.
    
    Parameters:
    -----------
    parent : mainWindow instance
        The parent window containing all display settings
    line : dict
        Dictionary containing line measurement information
        
    Returns:
    --------
    None
        Updates the 'text' key in the line dictionary
    """
    try:
        # Format text for equivalent width measurement
        EW_det_text = f"{line['EW']:.0f} \u00B1 {line['EWsig']:.0f} m\u00C5"
        
        # Format text for upper limit
        EW_limit_text = f"<{2.*line['EWsig']:.0f} m\u00C5"
        
        # Format text for column density
        if line['N'] > 0:
            log_N = np.log10(line['N'])
            log_N_err = np.log10(line['N'] + line['Nsig']) - log_N
            logN_det_text = f"{log_N:.2f} \u00B1 {log_N_err:.3f} /cm\u00B2"
        else:
            logN_det_text = "N/A"

        # Determine which text to display based on flags
        # line['flag'] is the line specific upper/lower/detections
        # parent.pFlag is the toggle button for which to show
        
        # Measurement: 
        if line['flag'] == 0:
            if parent.pFlag == 0:
                text = ""
            elif parent.pFlag == 1:
                text = EW_det_text
            elif parent.pFlag == 2:
                text = logN_det_text

        # Upper Limit
        elif line['flag'] == 1:
            if parent.pFlag == 0:
                text = ""
            elif parent.pFlag == 1:
                text = EW_limit_text
            elif parent.pFlag == 2:
                text = f"<{np.log10(line['Nsig']):.2f}"

        # Lower Limit
        elif line['flag'] == 2:
            if parent.pFlag == 0:
                text = ""
            elif parent.pFlag == 1:
                text = EW_det_text
            elif parent.pFlag == 2:
                text = f">{np.log10(line['N']):.2f}"
                
        line['text'] = text
        
    except Exception as e:
        print(f"Error generating display text: {e}")
        line['text'] = "Error"
