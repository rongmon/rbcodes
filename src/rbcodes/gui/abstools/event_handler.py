"""
Event handler module for the absorption line analysis toolbox.
This module handles all user interactions including mouse clicks,
key presses, and mouse motion events.
"""

import numpy as np
from PyQt5.QtWidgets import QInputDialog
from plotting import Plotting
from equivalent_width import EquivalentWidth

class EventHandler:
    """
    Class for handling all user interaction events in the application.
    """
    
    @staticmethod
    def on_motion(parent, event):
        """
        Handle mouse motion events to track which axes is active and
        provide visual feedback.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all data and UI elements
        event : matplotlib.backend_bases.MouseEvent
            The mouse motion event
            
        Returns:
        --------
        None
        """
        try:
            if event.xdata is not None and event.ydata is not None:
                # Determine the page and axis that contains the mouse pointer
                found = False
                for p in range(len(parent.axesL)):
                    for i in range(len(parent.axesL[p])):
                        if parent.axesL[p][i] is not None and event.inaxes == parent.axesL[p][i]:
                            parent.page = p
                            parent.Lidx = i
                            parent.Ridx = None
                            found = True
                            break
                        elif parent.axesR[p][i] is not None and event.inaxes == parent.axesR[p][i]:
                            parent.page = p
                            parent.Ridx = i
                            parent.Lidx = None
                            found = True
                            break
                    if found:
                        break
                
                # If we found an axis, update the highlighting
                if found:
                    # Determine which axis to highlight
                    if parent.Lidx is not None:
                        current_axis = parent.axesL[parent.page][parent.Lidx]
                    else:
                        current_axis = parent.axesR[parent.page][parent.Ridx]
                    
                    # If this is a different axis than was previously highlighted
                    if parent.old_axes and parent.old_axes != current_axis:
                        # Remove highlighting from old axis
                        for pos in ['top', 'bottom', 'left', 'right']:
                            parent.old_axes.spines[pos].set_edgecolor('black')
                            parent.old_axes.spines[pos].set_linewidth(0.5)
                        
                    # If this is a new axis to highlight
                    if parent.old_axes != current_axis:
                        # Add highlighting to new axis
                        for pos in ['top', 'bottom', 'left', 'right']:
                            current_axis.spines[pos].set_edgecolor('#01DF01')
                            current_axis.spines[pos].set_linewidth(2)
                        
                        # Update the reference to the currently highlighted axis
                        parent.old_axes = current_axis
                        
                        # Redraw the figure to show the highlighting
                        try:
                            parent.figs[parent.page].canvas.draw_idle()
                        except:
                            parent.figs[parent.page].canvas.draw()
                            
        except Exception as e:
            print(f"Error in motion event handler: {e}")
            
    @staticmethod
    def on_press(parent, event):
        """
        Handle key press events for various functionality controls.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all data and UI elements
        event : matplotlib.backend_bases.KeyEvent
            The key press event
            
        Returns:
        --------
        None
        """
        print(f"Key press detected: {event.key}")  # Debug output

        try:
            # Handle 'v' key for manual entry of limits
            if event.key == 'v':
                # Check if we're in the left panel
                if parent.Lidx is not None and parent.old_axes in parent.axesL[parent.page]:
                    # Manual mask limits for continuum fitting
                    mask_reg, ok = QInputDialog.getText(
                        parent, 'Manual Mask Limits', 
                        'Input Region to Mask (e.g. 200,250)'
                    )
                    
                    if ok:
                        key_idx = (parent.page * 6) + parent.Lidx
                        vel = parent.ions[parent.keys[key_idx]]['vel']
                        wc = parent.ions[parent.keys[key_idx]]['wc']
                        
                        try:
                            mask = mask_reg.split(',')
                            mask = np.array(mask).astype('float32')
                            
                            wc = ((vel < mask[0]) | (vel > mask[1])) & wc
                            parent.ions[parent.keys[key_idx]]['wc'] = wc
                            Plotting.plot(parent, parent.Lidx, modify=True)
                        except Exception as e:
                            print(f"Error processing mask limits: {e}")
                            
                # Check if we're in the right panel
                elif parent.Ridx is not None and parent.old_axes in parent.axesR[parent.page]:
                    # Manual EW integration limits
                    integ_lims, ok = QInputDialog.getText(
                        parent, 'Manual EW Limits', 
                        'Input integration region (eg -100,100)'
                    )
                    
                    if ok:
                        key_idx = (parent.page * 6) + parent.Ridx
                        
                        try:
                            integ_lims = integ_lims.split(',')
                            integ_lims = np.array(integ_lims).astype('float32')
                            
                            parent.ions[parent.keys[key_idx]]['EWlims'][0] = integ_lims[0]
                            parent.ions[parent.keys[key_idx]]['EWlims'][1] = integ_lims[1]
                            Plotting.plot(parent, parent.Ridx, modify=False, Print=False)
                        except Exception as e:
                            print(f"Error processing integration limits: {e}")
            
            # Handle 'up' arrow key for increasing polynomial order                
            elif event.key == 'up':
                if parent.Lidx is not None:
                    key_idx = parent.page * 6 + parent.Lidx
                    parent.ions[parent.keys[key_idx]]['order'] += 1
                    Plotting.plot(parent, parent.Lidx, modify=True)
                else:
                    print('Click on a left transition window first')
            
            # Handle 'down' arrow key for decreasing polynomial order
            elif event.key == 'down':
                if parent.Lidx is not None:
                    key_idx = parent.page * 6 + parent.Lidx
                    parent.ions[parent.keys[key_idx]]['order'] -= 1
                    if parent.ions[parent.keys[key_idx]]['order'] < 0:
                        parent.ions[parent.keys[key_idx]]['order'] = 0
                    Plotting.plot(parent, parent.Lidx, modify=True)
                else:
                    print('Click on a transition window first')
            
            # Handle 'm' key for measuring equivalent width of current subplot
            elif event.key == 'm':
                if parent.Ridx is not None:
                    key_idx = parent.page * 6 + parent.Ridx
                    EquivalentWidth.calculate(
                        parent, parent.page, parent.Ridx,
                        parent.ions[parent.keys[key_idx]]['EWlims']
                    )
                else:
                    print('Click on a right transition window first')
            
            # Handle 'M' key for measuring equivalent width of all subplots
            elif event.key == 'M':
                for jj in range(len(parent.figs)):
                    orig_page = parent.page
                    parent.page = jj
                    for ii in range(parent.nions[parent.page]):
                        EquivalentWidth.calculate(
                            parent, parent.page, ii,
                            parent.ions[parent.keys[ii + parent.page * 6]]['EWlims']
                        )
                    parent.page = orig_page
            
            # Handle 'V' key for applying current velocity limits to all ions
            elif event.key == 'V':
                if parent.Ridx is not None:
                    EWlims = parent.ions[parent.keys[parent.Ridx + 6 * parent.page]]['EWlims']
                    for jj in range(len(parent.figs)):
                        orig_page = parent.page
                        parent.page = jj
                        for ii in range(parent.nions[parent.page]):
                            parent.ions[parent.keys[ii + parent.page * 6]]['EWlims'] = EWlims
                            Plotting.plot(parent, ii, modify=False, Print=False)
                        parent.page = orig_page
                else:
                    print('Click on a right transition window first')
            
            # Handle numeric keys for flagging absorber detection type
            elif event.key in ['0', '1', '2']:
                if parent.Ridx is not None:
                    key_idx = parent.page * 6 + parent.Ridx
                    parent.ions[parent.keys[key_idx]]['flag'] = int(event.key)
                    Plotting.plot(parent, parent.Ridx, modify=False, Print=True)
                else:
                    print('Click on a right transition window first')
            
            # Handle 't' key for toggling EW/N display mode
            elif event.key == 't':
                parent.pFlag = (parent.pFlag + 1) % 3
                for jj in range(len(parent.figs)):
                    orig_page = parent.page
                    parent.page = jj
                    for ii in range(parent.nions[parent.page]):
                        Plotting.plot(parent, ii, modify=False, Print=True)
                    parent.page = orig_page
                    
            # Handle 'q' key for graceful exit
            elif event.key == 'q':
                print("Exit requested via 'q' key press")
                # Use the clean_exit method for a proper shutdown
                if hasattr(parent, 'clean_exit'):
                    parent.clean_exit()
                else:
                    # Fallback to regular close if clean_exit is not available
                    parent.close()
                    
        except Exception as e:
            print(f"Error in key press event handler: {e}")
            import traceback
            traceback.print_exc()  # Full stack trace for debugging
            
    @staticmethod
    def on_click(parent, event):
        """
        Handle mouse click events for selecting velocity regions and limits.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all data and UI elements
        event : matplotlib.backend_bases.MouseEvent
            The mouse click event
            
        Returns:
        --------
        None
        """
        try:
            if event.button in [1, 3]:  # Left or right click
                # Left hand side is for continuum fitting
                if parent.Lidx is not None:
                    key_idx = (parent.page * 6) + parent.Lidx
                    vel = parent.ions[parent.keys[key_idx]]['vel']
                    name = parent.ions[parent.keys[key_idx]]['name']
                    
                    # Handle first click to set first velocity limit
                    if parent.vclim is None:
                        parent.vclim = [event.xdata]
                        parent.axesL[parent.page][parent.Lidx].plot(event.xdata, event.ydata, 'ro', ms=5)
                        parent.figs[parent.page].canvas.draw()
                    else:
                        # Handle second click to set second velocity limit and update mask
                        vclim = np.sort(np.append(parent.vclim, event.xdata))
                        parent.vclim = None
                        wc = parent.ions[parent.keys[key_idx]]['wc']
                        
                        if event.button == 1:  # Left click adds to masks
                            wc = ((vel < vclim[0]) | (vel > vclim[1])) & wc
                        else:  # Right click removes mask
                            wc = ((vel > vclim[0]) & (vel < vclim[1])) | wc
                            
                        # Update wc for plotting
                        parent.ions[parent.keys[key_idx]]['wc'] = wc
                        
                        # Replot
                        Plotting.plot(parent, parent.Lidx, modify=True)
                
                # Right hand side for picking velocity limits for EW measurements
                if parent.Ridx is not None:
                    key_idx = (parent.page * 6) + parent.Ridx
                    
                    # If left click then define leftward vel limit
                    if event.button == 1:
                        parent.EWlim[0] = event.xdata  # Used for plotting all with same range 'V' command
                        parent.ions[parent.keys[key_idx]]['EWlims'][0] = event.xdata
                        # Plot selected limits
                        Plotting.plot(parent, parent.Ridx, modify=False, Print=False)
                        
                    # If right click define rightward vel limit
                    elif event.button == 3:
                        parent.EWlim[1] = event.xdata
                        parent.ions[parent.keys[key_idx]]['EWlims'][1] = event.xdata
                        Plotting.plot(parent, parent.Ridx, modify=False, Print=False)
                        
        except Exception as e:
            print(f"Error in mouse click event handler: {e}")
            parent.vclim = None  # Reset in case of error