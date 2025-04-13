"""
Resource cleanup module for the absorption line analysis toolbox.
This module provides functions for safely cleaning up resources and preventing
segmentation faults when closing the application.
"""

import sys
import gc
import weakref
from PyQt5.QtWidgets import QApplication

class ResourceCleanup:
    """
    Class for safely cleaning up resources to prevent segmentation faults.
    """
    
    @staticmethod
    def disconnect_matplotlib_events(window):
        """
        Safely disconnect all matplotlib event handlers.
        
        Parameters:
        -----------
        window : MainWindow instance
            The main window containing matplotlib connections
        """
        try:
            # Get all attributes that may be connection IDs
            for attr_name in dir(window):
                if attr_name.startswith('cid') and isinstance(getattr(window, attr_name), int):
                    cid = getattr(window, attr_name)
                    
                    # Try to disconnect from all canvases
                    for canvas in window.canvas:
                        try:
                            if hasattr(canvas, 'mpl_disconnect'):
                                canvas.mpl_disconnect(cid)
                                print(f"Disconnected event handler {attr_name}")
                        except Exception as e:
                            print(f"Error disconnecting {attr_name}: {e}")
        except Exception as e:
            print(f"Error in disconnect_matplotlib_events: {e}")
    
    @staticmethod
    def clear_matplotlib_figures(window):
        """
        Safely clear all matplotlib figures.
        
        Parameters:
        -----------
        window : MainWindow instance
            The main window containing matplotlib figures
        """
        try:
            import matplotlib.pyplot as plt
            
            # Clear and close all figures
            for i, fig in enumerate(window.figs):
                try:
                    # Remove all axes and clear the figure
                    if hasattr(fig, 'clear'):
                        fig.clear()
                    
                    # Close the figure
                    plt.close(fig)
                    print(f"Cleared figure {i}")
                except Exception as e:
                    print(f"Error clearing figure {i}: {e}")
        except Exception as e:
            print(f"Error in clear_matplotlib_figures: {e}")
    
    @staticmethod
    def disconnect_qt_signals(window):
        """
        Safely disconnect all Qt signals.
        
        Parameters:
        -----------
        window : MainWindow instance
            The main window containing Qt signals
        """
        try:
            # Try to disconnect MainWindowSignals
            if hasattr(window, 'signals'):
                # The following uses a try/except for each signal to ensure
                # one failure doesn't prevent disconnecting other signals
                
                # Note: The proper way to disconnect signals varies depending on your PyQt version
                # and how the connections were established. This is a general approach.
                
                signals = window.signals
                
                try:
                    if hasattr(signals, 'error_occurred'):
                        signals.error_occurred.disconnect()
                except Exception as e:
                    print(f"Error disconnecting error_occurred: {e}")
                
                try:
                    if hasattr(signals, 'status_message'):
                        signals.status_message.disconnect()
                except Exception as e:
                    print(f"Error disconnecting status_message: {e}")
                
                try:
                    if hasattr(signals, 'tab_changed'):
                        signals.tab_changed.disconnect()
                except Exception as e:
                    print(f"Error disconnecting tab_changed: {e}")
                
                try:
                    if hasattr(signals, 'data_updated'):
                        signals.data_updated.disconnect()
                except Exception as e:
                    print(f"Error disconnecting data_updated: {e}")
                
                try:
                    if hasattr(signals, 'file_loaded'):
                        signals.file_loaded.disconnect()
                except Exception as e:
                    print(f"Error disconnecting file_loaded: {e}")
                
                try:
                    if hasattr(signals, 'file_saved'):
                        signals.file_saved.disconnect()
                except Exception as e:
                    print(f"Error disconnecting file_saved: {e}")
        except Exception as e:
            print(f"Error in disconnect_qt_signals: {e}")
    

    @staticmethod
    def safe_exit(window):
        """
        Perform a safe exit of the application to prevent segmentation faults.
        
        Parameters:
        -----------
        window : MainWindow instance
            The main window to properly close
        """
        try:
            print("Starting safe application exit...")
            
            # Step 1: Disconnect all matplotlib events
            ResourceCleanup.disconnect_matplotlib_events(window)
            
            # Step 2: Clear and close matplotlib figures
            ResourceCleanup.clear_matplotlib_figures(window)
            
            # Step 3: Disconnect all Qt signals
            ResourceCleanup.disconnect_qt_signals(window)
            
            # Step 4: Clear references to large objects
            try:
                for attr_name in ['axesL', 'axesR', 'canvas', 'figs', 'ions']:
                    if hasattr(window, attr_name):
                        setattr(window, attr_name, None)
                print("Cleared large object references")
            except Exception as e:
                print(f"Error clearing references: {e}")
            
            # Step 5: Force a garbage collection run
            gc.collect()
            print("Garbage collection complete")
            
            # Step 6: IMPORTANT - Use os._exit instead of sys.exit
            # This bypasses all Python cleanup which might be triggering the segfault
            print("Exiting application immediately...")
            import os
            os._exit(0)  # Force immediate termination with no cleanup
                    
        except Exception as e:
            print(f"Error during safe exit: {e}")
            # If all else fails, force an immediate exit
            import os
            os._exit(1)    