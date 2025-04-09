"""
Cleanup module for the absorption line analysis toolbox.
This module handles proper resource cleanup to prevent memory leaks
and segmentation faults on application exit.
"""

import gc
import sys
import weakref
import numpy as np

class ResourceCleanup:
    """
    Class for managing resource cleanup and preventing segmentation faults.
    """
    
    @staticmethod
    def disconnect_mpl_callbacks(figure):
        """
        Disconnect all matplotlib callbacks from a figure safely.
        
        Parameters:
        -----------
        figure : matplotlib.figure.Figure
            The figure to disconnect callbacks from
            
        Returns:
        --------
        None
        """
        if hasattr(figure, 'canvas') and figure.canvas is not None:
            # Store all callback IDs
            callback_ids = []
            
            # Get callbacks from the canvas
            if hasattr(figure.canvas, 'callbacks') and figure.canvas.callbacks is not None:
                try:
                    for cid in figure.canvas.callbacks.callbacks:
                        callback_ids.append(cid)
                except:
                    # If we can't get the callbacks, just move on
                    pass
                    
            # Disconnect each callback
            for cid in callback_ids:
                try:
                    figure.canvas.mpl_disconnect(cid)
                except:
                    # Ignore errors when disconnecting
                    pass
    
    @staticmethod
    def cleanup_figures(figures):
        """
        Clean up matplotlib figures safely without causing segmentation faults.
        
        Parameters:
        -----------
        figures : list
            List of matplotlib figures to clean up
            
        Returns:
        --------
        None
        """
        if figures is not None:
            for fig in figures:
                if fig is not None:
                    # Disconnect all callbacks first
                    ResourceCleanup.disconnect_mpl_callbacks(fig)
                    
                    # Close the figure directly without using pyplot
                    try:
                        # Clear the figure
                        fig.clear()
                        
                        # Close the canvas if it exists
                        if hasattr(fig, 'canvas') and fig.canvas is not None:
                            fig.canvas.close()
                    except Exception as e:
                        print(f"Warning during figure cleanup: {e}")
            
            # Clear the figures list
            figures.clear()
    
    @staticmethod
    def cleanup_qt_widgets(widgets):
        """
        Clean up PyQt widgets safely.
        
        Parameters:
        -----------
        widgets : list
            List of PyQt widgets to clean up
            
        Returns:
        --------
        None
        """
        if widgets is not None:
            for widget in widgets:
                if widget is not None:
                    try:
                        # Schedule for deletion later
                        widget.deleteLater()
                    except:
                        # Ignore errors in widget cleanup
                        pass
            
            # Clear the widgets list
            widgets.clear()
    
    @staticmethod
    def cleanup_application(parent):
        """
        Clean up all resources associated with the application.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all resources
            
        Returns:
        --------
        None
        """
        try:
            # Disconnect all matplotlib callbacks first
            if hasattr(parent, 'cid1'):
                for attr_name in dir(parent):
                    if attr_name.startswith('cid') and hasattr(parent, attr_name):
                        cid = getattr(parent, attr_name)
                        if isinstance(cid, int):
                            for fig in parent.figs:
                                try:
                                    fig.canvas.mpl_disconnect(cid)
                                except:
                                    pass
            
            # Clean up large data arrays to free memory
            if hasattr(parent, 'keys') and hasattr(parent, 'ions'):
                for key in parent.keys:
                    if key in parent.ions:
                        for item in ['vel', 'wave', 'flux', 'error', 'weight', 'cont']:
                            if item in parent.ions[key] and isinstance(parent.ions[key][item], np.ndarray):
                                parent.ions[key][item] = np.array([])
            
            # Clean up matplotlib figures
            if hasattr(parent, 'figs'):
                ResourceCleanup.cleanup_figures(parent.figs)
            
            # Clean up Qt widgets
            if hasattr(parent, 'tabs'):
                ResourceCleanup.cleanup_qt_widgets(parent.tabs)
            
            # Clean up canvas objects
            if hasattr(parent, 'canvas'):
                ResourceCleanup.cleanup_qt_widgets(parent.canvas)
            
            # Clean up axes
            if hasattr(parent, 'axesL'):
                parent.axesL = None
            if hasattr(parent, 'axesR'):
                parent.axesR = None
            parent.old_axes = None
            
            # Force garbage collection
            gc.collect()
            
        except Exception as e:
            print(f"Error during cleanup: {e}")

    @staticmethod
    def register_cleanup_on_exit(app, parent):
        """
        Register cleanup functions to be called on application exit.
        
        Parameters:
        -----------
        app : QApplication
            The Qt application instance
        parent : mainWindow instance
            The parent window containing all resources
            
        Returns:
        --------
        None
        """
        from PyQt5.QtCore import QObject, pyqtSignal, QTimer
        
        # Create a cleanup handler
        class CleanupHandler(QObject):
            aboutToQuit = pyqtSignal()
            
            def __init__(self, parent):
                super().__init__()
                self.parent_ref = weakref.ref(parent)
                
            def cleanup(self):
                parent = self.parent_ref()
                if parent is not None:
                    try:
                        ResourceCleanup.cleanup_application(parent)
                        
                        # Schedule application termination after a brief delay
                        QTimer.singleShot(200, lambda: sys.exit(0))
                    except Exception as e:
                        print(f"Error during cleanup: {e}")
                        sys.exit(0)  # Exit anyway
        
        # Create and connect the handler
        handler = CleanupHandler(parent)
        app.aboutToQuit.connect(handler.cleanup)
        
        # Keep a reference to the handler to prevent it from being garbage collected
        parent._cleanup_handler = handler