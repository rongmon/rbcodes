"""
Event handler module for the absorption line analysis toolbox.
This module handles all user interactions including mouse clicks,
key presses, and mouse motion events.
Improved with signal-slot communication and better error handling.
"""

import numpy as np
from PyQt5.QtWidgets import QInputDialog
from rbcodes.GUIs.abstools.plotting import Plotting
from rbcodes.GUIs.abstools.equivalent_width import EquivalentWidth

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
            # We don't signal errors for motion events since they happen frequently
            
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
                            
                            # Signal that mask has been updated
                            if hasattr(parent, 'signals'):
                                parent.signals.data_updated.emit({
                                    "action": "mask_updated",
                                    "ion_index": key_idx,
                                    "mask_limits": mask.tolist()
                                })
                            
                            Plotting.plot(parent, parent.Lidx, modify=True)
                            
                            if hasattr(parent, 'signals'):
                                parent.signals.status_message.emit(f"Mask region set: {mask[0]}-{mask[1]}")
                                
                        except Exception as e:
                            print(f"Error processing mask limits: {e}")
                            if hasattr(parent, 'signals'):
                                parent.signals.error_occurred.emit(f"Invalid mask format: {str(e)}")
                            
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
                            
                            # Signal that EW limits have been updated
                            if hasattr(parent, 'signals'):
                                parent.signals.data_updated.emit({
                                    "action": "ew_limits_updated",
                                    "ion_index": key_idx,
                                    "ew_limits": integ_lims.tolist()
                                })
                                
                            Plotting.plot(parent, parent.Ridx, modify=False, Print=False)
                            
                            if hasattr(parent, 'signals'):
                                parent.signals.status_message.emit(f"Integration limits set: {integ_lims[0]}-{integ_lims[1]}")
                                
                        except Exception as e:
                            print(f"Error processing integration limits: {e}")
                            if hasattr(parent, 'signals'):
                                parent.signals.error_occurred.emit(f"Invalid integration limit format: {str(e)}")
            
            # Handle 'up' arrow key for increasing polynomial order                
            elif event.key == 'up':
                if parent.Lidx is not None:
                    key_idx = parent.page * 6 + parent.Lidx
                    parent.ions[parent.keys[key_idx]]['order'] += 1
                    
                    # Signal that polynomial order has changed
                    if hasattr(parent, 'signals'):
                        parent.signals.data_updated.emit({
                            "action": "poly_order_changed",
                            "ion_index": key_idx,
                            "new_order": parent.ions[parent.keys[key_idx]]['order']
                        })
                    
                    Plotting.plot(parent, parent.Lidx, modify=True)
                    
                    if hasattr(parent, 'signals'):
                        parent.signals.status_message.emit(f"Polynomial order increased to {parent.ions[parent.keys[key_idx]]['order']}")
                else:
                    if hasattr(parent, 'signals'):
                        parent.signals.status_message.emit('Click on a left transition window first')
                    else:
                        print('Click on a left transition window first')
            
            # Handle 'down' arrow key for decreasing polynomial order
            elif event.key == 'down':
                if parent.Lidx is not None:
                    key_idx = parent.page * 6 + parent.Lidx
                    parent.ions[parent.keys[key_idx]]['order'] -= 1
                    if parent.ions[parent.keys[key_idx]]['order'] < 0:
                        parent.ions[parent.keys[key_idx]]['order'] = 0
                    
                    # Signal that polynomial order has changed
                    if hasattr(parent, 'signals'):
                        parent.signals.data_updated.emit({
                            "action": "poly_order_changed",
                            "ion_index": key_idx,
                            "new_order": parent.ions[parent.keys[key_idx]]['order']
                        })
                    
                    Plotting.plot(parent, parent.Lidx, modify=True)
                    
                    if hasattr(parent, 'signals'):
                        parent.signals.status_message.emit(f"Polynomial order decreased to {parent.ions[parent.keys[key_idx]]['order']}")
                else:
                    if hasattr(parent, 'signals'):
                        parent.signals.status_message.emit('Click on a transition window first')
                    else:
                        print('Click on a transition window first')
            
            # Handle 'm' key for measuring equivalent width of current subplot
            elif event.key == 'm':
                if parent.Ridx is not None:
                    key_idx = parent.page * 6 + parent.Ridx
                    success = EquivalentWidth.calculate(
                        parent, parent.page, parent.Ridx,
                        parent.ions[parent.keys[key_idx]]['EWlims']
                    )
                    
                    if success and hasattr(parent, 'signals'):
                        # Get the calculated values
                        ion_name = parent.keys[key_idx]
                        ew_value = parent.ions[ion_name]['EW']
                        ew_error = parent.ions[ion_name]['EWsig']
                        
                        # Signal that EW has been measured
                        parent.signals.ew_measured.emit(ion_name, ew_value, ew_error)
                        parent.signals.status_message.emit(f"Measured EW for {ion_name}: {ew_value:.2f} ± {ew_error:.2f} mÅ")
                else:
                    if hasattr(parent, 'signals'):
                        parent.signals.status_message.emit('Click on a right transition window first')
                    else:
                        print('Click on a right transition window first')
            
            # Handle 'M' key for measuring equivalent width of all subplots
            elif event.key == 'M':
                measured_count = 0
                total_ions = 0
                
                for jj in range(len(parent.figs)):
                    orig_page = parent.page
                    parent.page = jj
                    for ii in range(parent.nions[parent.page]):
                        total_ions += 1
                        key_idx = ii + parent.page * 6
                        success = EquivalentWidth.calculate(
                            parent, parent.page, ii,
                            parent.ions[parent.keys[key_idx]]['EWlims']
                        )
                        if success:
                            measured_count += 1
                            
                            if hasattr(parent, 'signals'):
                                # Get the calculated values
                                ion_name = parent.keys[key_idx]
                                ew_value = parent.ions[ion_name]['EW']
                                ew_error = parent.ions[ion_name]['EWsig']
                                
                                # Signal that EW has been measured
                                parent.signals.ew_measured.emit(ion_name, ew_value, ew_error)
                        
                    parent.page = orig_page
                
                if hasattr(parent, 'signals'):
                    if measured_count == total_ions:
                        parent.signals.status_message.emit(f"Measured EW for all {measured_count} transitions")
                    else:
                        parent.signals.status_message.emit(f"Measured EW for {measured_count} of {total_ions} transitions")
            
            # Handle 'V' key for applying current velocity limits to all ions
            elif event.key == 'V':
                if parent.Ridx is not None:
                    key_idx = parent.Ridx + 6 * parent.page
                    EWlims = parent.ions[parent.keys[key_idx]]['EWlims']
                    
                    # Signal that EW limits will be applied to all ions
                    if hasattr(parent, 'signals'):
                        parent.signals.data_updated.emit({
                            "action": "ew_limits_global_update",
                            "ew_limits": EWlims
                        })
                    
                    for jj in range(len(parent.figs)):
                        orig_page = parent.page
                        parent.page = jj
                        for ii in range(parent.nions[parent.page]):
                            parent.ions[parent.keys[ii + parent.page * 6]]['EWlims'] = EWlims
                            Plotting.plot(parent, ii, modify=False, Print=False)
                        parent.page = orig_page
                    
                    if hasattr(parent, 'signals'):
                        parent.signals.status_message.emit(f"Applied velocity limits {EWlims[0]:.1f}-{EWlims[1]:.1f} to all transitions")
                else:
                    if hasattr(parent, 'signals'):
                        parent.signals.status_message.emit('Click on a right transition window first')
                    else:
                        print('Click on a right transition window first')
            
            # Handle numeric keys for flagging absorber detection type
            elif event.key in ['0', '1', '2']:
                if parent.Ridx is not None:
                    key_idx = parent.page * 6 + parent.Ridx
                    flag_value = int(event.key)
                    parent.ions[parent.keys[key_idx]]['flag'] = flag_value
                    
                    # Signal that the flag has been changed
                    if hasattr(parent, 'signals'):
                        flag_types = {0: "Detection", 1: "Upper Limit", 2: "Lower Limit"}
                        parent.signals.data_updated.emit({
                            "action": "flag_updated",
                            "ion_index": key_idx,
                            "flag": flag_value,
                            "flag_type": flag_types.get(flag_value, "Unknown")
                        })
                        parent.signals.status_message.emit(f"Set {parent.keys[key_idx]} as {flag_types.get(flag_value, 'Unknown')}")
                    
                    Plotting.plot(parent, parent.Ridx, modify=False, Print=True)
                else:
                    if hasattr(parent, 'signals'):
                        parent.signals.status_message.emit('Click on a right transition window first')
                    else:
                        print('Click on a right transition window first')
            
            # Handle 't' key for toggling EW/N display mode
            elif event.key == 't':
                parent.pFlag = (parent.pFlag + 1) % 3
                display_modes = ["None", "Equivalent Width", "Column Density"]
                
                if hasattr(parent, 'signals'):
                    parent.signals.data_updated.emit({
                        "action": "display_mode_changed",
                        "display_mode": parent.pFlag,
                        "display_name": display_modes[parent.pFlag]
                    })
                    parent.signals.status_message.emit(f"Display mode changed to: {display_modes[parent.pFlag]}")
                
                for jj in range(len(parent.figs)):
                    orig_page = parent.page
                    parent.page = jj
                    for ii in range(parent.nions[parent.page]):
                        Plotting.plot(parent, ii, modify=False, Print=True)
                    parent.page = orig_page
                    
            # Handle 'q' key for graceful exit
            elif event.key == 'q':
                print("Exit requested via 'q' key press")
                if hasattr(parent, 'signals'):
                    parent.signals.status_message.emit("Exiting application...")
                
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
            
            if hasattr(parent, 'signals'):
                parent.signals.error_occurred.emit(f"Error processing keyboard command: {str(e)}")
            
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
                        
                        if hasattr(parent, 'signals'):
                            parent.signals.status_message.emit(f"Set first continuum mask point at {event.xdata:.1f}")
                    else:
                        # Handle second click to set second velocity limit and update mask
                        vclim = np.sort(np.append(parent.vclim, event.xdata))
                        parent.vclim = None
                        wc = parent.ions[parent.keys[key_idx]]['wc']
                        
                        if event.button == 1:  # Left click adds to masks
                            wc = ((vel < vclim[0]) | (vel > vclim[1])) & wc
                            action_type = "added"
                        else:  # Right click removes mask
                            wc = ((vel > vclim[0]) & (vel < vclim[1])) | wc
                            action_type = "removed"
                            
                        # Update wc for plotting
                        parent.ions[parent.keys[key_idx]]['wc'] = wc
                        
                        # Signal that the mask has been updated
                        if hasattr(parent, 'signals'):
                            parent.signals.data_updated.emit({
                                "action": "mask_updated",
                                "ion_index": key_idx,
                                "mask_limits": vclim.tolist(),
                                "mask_action": action_type
                            })
                            parent.signals.status_message.emit(f"{action_type.capitalize()} mask region: {vclim[0]:.1f}-{vclim[1]:.1f}")
                        
                        # Replot
                        Plotting.plot(parent, parent.Lidx, modify=True)
                
                # Right hand side for picking velocity limits for EW measurements
                if parent.Ridx is not None:
                    key_idx = (parent.page * 6) + parent.Ridx
                    
                    # If left click then define leftward vel limit
                    if event.button == 1:
                        parent.EWlim[0] = event.xdata  # Used for plotting all with same range 'V' command
                        parent.ions[parent.keys[key_idx]]['EWlims'][0] = event.xdata
                        
                        # Signal that EW limit has been updated
                        if hasattr(parent, 'signals'):
                            parent.signals.data_updated.emit({
                                "action": "ew_limit_updated",
                                "ion_index": key_idx,
                                "limit_type": "left",
                                "value": event.xdata
                            })
                            parent.signals.status_message.emit(f"Set left integration limit to {event.xdata:.1f}")
                        
                        # Plot selected limits
                        Plotting.plot(parent, parent.Ridx, modify=False, Print=False)
                        
                    # If right click define rightward vel limit
                    elif event.button == 3:
                        parent.EWlim[1] = event.xdata
                        parent.ions[parent.keys[key_idx]]['EWlims'][1] = event.xdata
                        
                        # Signal that EW limit has been updated
                        if hasattr(parent, 'signals'):
                            parent.signals.data_updated.emit({
                                "action": "ew_limit_updated",
                                "ion_index": key_idx,
                                "limit_type": "right",
                                "value": event.xdata
                            })
                            parent.signals.status_message.emit(f"Set right integration limit to {event.xdata:.1f}")
                        
                        Plotting.plot(parent, parent.Ridx, modify=False, Print=False)
                        
        except Exception as e:
            print(f"Error in mouse click event handler: {e}")
            if hasattr(parent, 'signals'):
                parent.signals.error_occurred.emit(f"Error processing mouse click: {str(e)}")
            parent.vclim = None  # Reset in case of error