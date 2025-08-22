import numpy as np
import pandas as pd
import traceback
import sys

# Import required GUI libraries
from PyQt5.QtWidgets import QInputDialog, QDialog, QVBoxLayout
from PyQt5 import QtCore

# Import matplotlib libraries
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

# Import required modules with error handling
try:
    from rbcodes.IGM import rb_setline as line
    print("Imported rb_setline successfully")
except ImportError as e:
    print(f"Failed to import rb_setline: {e}")

try:
    from rbcodes.GUIs.abstools import Absorber as A
    print("Imported Absorber successfully")
except ImportError as e:
    print(f"Failed to import Absorber: {e}")

try:
    from rbcodes.utils import rb_utility as rt
    print("Imported rb_utility successfully")
    # Get color dictionary
    clr = rt.rb_set_color()
except ImportError as e:
    print(f"Failed to import rb_utility: {e}")
    # Fallback color dictionary if import fails
    clr = {
        'white': '#FFFFFF',
        'teal': '#008080',
        'red': '#FF0000',
        'yellow': '#FFFF00',
        'cyan': '#00FFFF',
        'orange2': '#FFA500',
        'pale_red': '#FF6666',
        'light_gray': '#CCCCCC'
    }

def prepare_absorber_object(z_abs, wave, flux, error, line_flg='LLS', vlim=[-1000, 1000]):
    """
    Prepare absorption line data for plotting
    
    Parameters:
    -----------
    z_abs : float
        Redshift of the absorber
    wave : array
        Wavelength array
    flux : array
        Flux array
    error : array
        Error array
    line_flg : str, optional
        Line list flag (default: 'LLS')
    vlim : list, optional
        Velocity limits (default: [-1000, 1000])
        
    Returns:
    --------
    dict
        Dictionary of ions data
    """
    # Read the full linelist
    data = line.read_line_list(line_flg)
    wavelist = []
    for i in range(0, len(data)):
        wavelist.append(data[i]['wrest'])
    
    wavelist = np.array(wavelist)
    
    # Select the lines within the wavelength range only
    q = np.where((wavelist > np.min(wave)/(1.+z_abs)) & (wavelist < (np.max(wave)/(1.+z_abs))))
    
    # Total transitions visible within the wavelength window
    nTot = len(q[0])
    
    wavelist_selected = wavelist[q]
    
    # Create absorber object
    absys = A.Absorber(z_abs, wave, flux, error, list(wavelist_selected), window_lim=vlim, nofrills=True)   
    
    return absys.ions

class vStack:
    def __init__(self, parent, wave, flux, error, line_flg, zabs=0, vlim=[-1000., 1000.]):
        """
        Initialize vStack class for displaying velocity plots of absorption lines
        
        Parameters:
        -----------
        parent : object
            Parent object (usually SpectralPlot)
        wave : array
            Wavelength array
        flux : array
            Flux array
        error : array
            Error array
        line_flg : str
            Line list flag
        zabs : float, optional
            Redshift of the absorber (default: 0)
        vlim : list, optional
            Velocity limits (default: [-1000., 1000.])
        """
        self.parent = parent
        self.zabs = zabs
        self.vlim = vlim
        
        print(f"Initializing vStack with z={zabs:.6f}, line_flg={line_flg}, vlim={vlim}")
        
        # Initialize ions data
        try:
            self.ions = prepare_absorber_object(zabs, wave, flux, error, line_flg=line_flg, vlim=vlim)
            print(f"prepare_absorber_object completed successfully")
        except Exception as e:
            print(f"Error in prepare_absorber_object: {str(e)}")
            traceback.print_exc()
            raise
        
        # Extract spectral properties
        self.z = self.ions['Target']['z'] 
        self.flux = self.ions['Target']['flux']
        self.wave = self.ions['Target']['wave'] 
        self.error = self.ions['Target']['error']
        
        # Get list of transitions (ions)
        self.keys = list(self.ions.keys())[:-1]  # last item is the full target spectrum
        self.nions = int(len(self.keys))
        
        print(f"Found {self.nions} transitions in wavelength range")
        
        # Set default flags for all ions (0 = non-detection)
        for i in self.keys:
            self.ions[i]['flag'] = 0
        
        # Configure display grid layout
        self.page = 1
        self.plotppage = 12  # plots per page
        self.nrow = int(self.plotppage/3)
        self.ncol = int((self.plotppage)/self.nrow)
        self.npages = int((self.nions//self.plotppage))
        
        # Calculate number of pages needed
        if self.nions % self.plotppage > 0:
            self.npages += 1
        
        print(f"Grid layout: {self.nrow}x{self.ncol}, {self.npages} pages")
        
        # Create the figure and axes
        self.fig = Figure(figsize=(12, 8))
        self.axes = list(range(self.plotppage))
        for i in range(self.plotppage):
            self.axes[i] = self.fig.add_subplot(self.nrow, self.ncol, i+1)
        self.axes = np.array(self.axes)
        
        # Initial plot
        self.vPlot()
        self.fig.subplots_adjust(hspace=0)
        
        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()
        
        # Connect keyboard events
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self.onkb)
        
        # Find the right parent widget and layout
        canvas_replaced = False
        
        try:
            # Get the real parent window from the canvas
            main_window = None
            
            # Try to find the MainWindow (parent_window)
            if hasattr(parent, 'parent_window'):
                main_window = parent.parent_window
                print("Found parent_window attribute on canvas")
            else:
                print("parent_window attribute not found on canvas, trying other approaches")
                
                # Check if parent itself is MainWindow
                from PyQt5.QtWidgets import QMainWindow
                if isinstance(parent, QMainWindow):
                    main_window = parent
                    print("Parent is QMainWindow")
            
            if main_window is not None:
                # Try different layout finding strategies
                if hasattr(main_window, 'right_layout'):
                    print("Found right_layout")
                    main_window.right_layout.removeWidget(parent)
                    main_window.right_layout.insertWidget(1, self.canvas)
                    canvas_replaced = True
                elif hasattr(main_window, 'main_splitter'):
                    print("Found main_splitter, looking for canvas in children")
                    
                    # Get the right widget in the splitter (usually index 1)
                    right_widget = main_window.main_splitter.widget(1)
                    
                    if hasattr(right_widget, 'layout'):
                        right_layout = right_widget.layout()
                        
                        # Find parent's canvas in the layout
                        for i in range(right_layout.count()):
                            item = right_layout.itemAt(i)
                            if item is not None and item.widget() == parent:
                                right_layout.removeWidget(parent)
                                right_layout.insertWidget(i, self.canvas)
                                canvas_replaced = True
                                print(f"Replaced canvas in right_widget layout at index {i}")
                                break
                            elif item is not None and hasattr(item.widget(), 'canvas') and item.widget().canvas == parent:
                                item_widget = item.widget()
                                right_layout.removeWidget(item_widget)
                                right_layout.insertWidget(i, self.canvas)
                                canvas_replaced = True
                                print(f"Replaced widget containing canvas at index {i}")
                                break
                else:
                    print("Could not find the right layout in MainWindow.")
                    print("MainWindow attributes:", [attr for attr in dir(main_window) if not attr.startswith('_')])
            else:
                print("Could not find MainWindow")
        except Exception as e:
            print(f"Error replacing canvas: {str(e)}")
            traceback.print_exc()
        
        # If we couldn't replace the canvas in the layout, show in a separate window
        if not canvas_replaced:
            print("Could not replace canvas in layout, showing in separate window")
            self.separate_window = QDialog()
            self.separate_window.setWindowTitle(f"vStack - z={zabs:.6f}")
            self.separate_window.setMinimumSize(800, 600)
            
            layout = QVBoxLayout(self.separate_window)
            layout.addWidget(self.canvas)
            
            self.separate_window.show()
        
        print("vStack initialization complete")

    def onkb(self, event):
        """
        Handle keyboard events
        
        Parameters:
        -----------
        event : matplotlib.backend_bases.KeyEvent
            Keyboard event
        """
        
        # Set up custom y-limit for an individual panel
        if event.key == 'Y':
            if event.inaxes in self.axes:
                i = np.where(event.inaxes == self.axes)[0][0] + self.plotppage*(self.page-1)
                Windowname = 'Manual y-Limits'
                instruction_text = 'Input range (e.g. 0.,2.)'
                ylim, ok = QInputDialog.getText(self.parent, Windowname, instruction_text)
                if ok:
                    try:
                        ylimit = ylim.split(',')
                        ylimit = np.array(ylimit).astype('float32')
                        self.vPlot(ploti=i, yrange=[ylimit[0], ylimit[1]])
                        print(f"Set y-limits to {ylimit} for panel {i}")
                    except Exception as e:
                        print(f"Error setting y-limits: {str(e)}")
        
        # Page right
        elif event.key == '>':
            self.page += 1
            if self.page > self.npages: 
                self.page = 1
            self.vPlot(clearpage=True)
            print(f"Moved to page {self.page} of {self.npages}")
        
        # Page left
        elif event.key == '<':
            self.page -= 1
            if self.page < 1: 
                self.page = self.npages
            self.vPlot(clearpage=True)
            print(f"Moved to page {self.page} of {self.npages}")
        
        # Toggle between detection-non-detection or blended    
        elif event.key == 'w':
            if event.inaxes in self.axes:
                i = np.where(event.inaxes == self.axes)[0][0] + self.plotppage*(self.page-1)
                # Toggle through flag states (0->1->2->3->0)
                temp_flag = self.ions[self.keys[i]]['flag']
                temp_flag = (temp_flag + 1) % 4
                self.ions[self.keys[i]]['flag'] = temp_flag
                self.vPlot(ploti=i, comment=False)
                print(f"Changed flag for {self.keys[i]} to {temp_flag} ({self.plotText(flag=temp_flag)})")
        
        # Save linelist and return to main display
        elif event.key == 'S':
            # Collect all lines with flags indicating detection (1, 2, or 3)
            lines_added = 0
            for key in self.keys:
                if self.ions[key]['flag'] == 0:  # Skip non-detections
                    continue
                
                wave_obs = self.ions[key]['lam_0_z']
                name = self.ions[key]['name']
                zabs = self.ions['Target']['z']
                
                # Format name based on flag
                display_name = name
                if self.ions[key]['flag'] == 2:
                    display_name += " [b]"  # Blended
                elif self.ions[key]['flag'] == 3:
                    display_name += " [p]"  # Potential/low confidence
                
                # Add to parent's line list
                new_row = pd.Series(data={'Name': display_name, 'Wave_obs': wave_obs, 'Zabs': zabs})
                
                # Check if parent has line_list attribute
                if hasattr(self.parent, 'line_list'):
                    self.parent.line_list = self.parent.line_list.append(new_row, ignore_index=True)
                    lines_added += 1
                elif hasattr(self.parent, 'canvas') and hasattr(self.parent.canvas, 'line_list'):
                    self.parent.canvas.line_list = self.parent.canvas.line_list.append(new_row, ignore_index=True)
                    lines_added += 1
                else:
                    print("Could not find line_list attribute on parent or parent.canvas")
            
            print(f"Added {lines_added} lines to line_list")
            
            # Close vStack and return to main display
            try:
                if hasattr(self, 'separate_window') and self.separate_window is not None:
                    self.separate_window.close()
                else:
                    # Find the original canvas
                    main_window = self.parent.parent_window
                    
                    # Try to restore original canvas
                    if hasattr(main_window, 'right_layout'):
                        main_window.right_layout.removeWidget(self.canvas)
                        main_window.right_layout.insertWidget(1, self.parent)
                    elif hasattr(main_window, 'main_splitter'):
                        right_widget = main_window.main_splitter.widget(1)
                        if hasattr(right_widget, 'layout'):
                            layout = right_widget.layout()
                            for i in range(layout.count()):
                                if layout.itemAt(i).widget() == self.canvas:
                                    layout.removeWidget(self.canvas)
                                    layout.insertWidget(i, self.parent)
                                    break
                
                # Close/delete the canvas
                self.canvas.close()
                del self.canvas
                
                print("Closed vStack and restored original canvas")
            except Exception as e:
                print(f"Error restoring original canvas: {str(e)}")
                traceback.print_exc()

    def vPlot(self, ploti=None, comment=False, yrange=None, clearpage=False):
        """
        Plot velocity panels
        
        Parameters:
        -----------
        ploti : int or array, optional
            Indices of panels to plot
        comment : str, optional
            Comment to add to plots
        yrange : list, optional
            Y-axis limits
        clearpage : bool, optional
            Whether to clear the page
        """
        if ploti is None:
            ploti = np.arange(self.plotppage*(self.page-1), min(self.plotppage*self.page, self.nions))
        else:
            ploti = [ploti]

        # Clear the axes if needed
        if clearpage:
            for i in range(self.plotppage):
                self.clearstuff(i)
                    
        # Plot the panels
        for i in ploti:
            if i < self.nions:  # Only plot if within range
                self.plotstuff(i, comment=comment, yrange=yrange)

        self.fig.canvas.draw()

    def clearstuff(self, i):
        """
        Clear a specific panel
        
        Parameters:
        -----------
        i : int
            Panel index
        """
        ax = self.axes[i % self.plotppage]
        ax.clear()
        ax.set_axis_off()

    def plotstuff(self, i, comment=False, yrange=False):
        """
        Plot a specific velocity panel
        
        Parameters:
        -----------
        i : int
            Panel index
        comment : str, optional
            Comment to add to plot
        yrange : list, optional
            Y-axis limits
        """
        ax = self.axes[i % self.plotppage]
        
        # Define variables for readability
        vel = self.ions[self.keys[i]]['vel']
        wave = self.ions[self.keys[i]]['wave']
        error = self.ions[self.keys[i]]['error']
        flux = self.ions[self.keys[i]]['flux']
        name = self.ions[self.keys[i]]['name']
        window_lim = self.ions[self.keys[i]]['window_lim']
        flag = self.ions[self.keys[i]]['flag']
        f0 = self.ions[self.keys[i]]['f']

        ax.clear()
        
        # Set a title giving the page number
        if i % self.plotppage == 0:
            ax.set_title('Page ' + str(self.page) + ' of ' + str(self.npages), color=clr['teal'])

        # Plot normalized flux and error
        ax.step(vel, flux/np.nanmean(flux), where='mid', color=clr['teal'],lw=0.5)
        ax.step(vel, error/np.nanmean(flux), where='mid', color=clr['orange2'],lw=0.5)
        
        # Add reference lines
        ax.axhline(1, color=clr['light_gray'], linestyle='dotted')
        ax.axvline(0, color=clr['light_gray'], linestyle='dotted')
        
        # Add transition name and oscillator strength
        ax.text(x=0.05, y=0.815, s=name, fontsize=8, transform=ax.transAxes, color=clr['red'])
        ax.text(x=0.75, y=0.815, s='f0: '+str(f0), fontsize=8, transform=ax.transAxes, color=clr['red'])
        
        # Add custom comment if provided
        if comment != False:
            ax.text(x=0.85, y=0.815, s=comment, fontsize=12, transform=ax.transAxes, color=clr['teal'])
            
        # Set custom y-range if provided
        if yrange != False:
            ax.set_ylim(yrange)
            
        # Set x-limits from velocity window
        ax.set_xlim(self.vlim)
        
        # Add detection status text
        if flag is not None:
            textout = self.plotText(flag=flag)
            
            if flag == 1:
                textcolor = clr['yellow']
            elif flag == 2:
                textcolor = clr['cyan']
            elif flag == 3:
                textcolor = clr['pale_red']
            else:
                textcolor = clr['light_gray']

            ax.text(x=0.05, y=0.01, s=textout, fontsize=12, transform=ax.transAxes, color=textcolor)

        # Calculate grid position
        panel_in_page = i % self.plotppage  # Position within current page (0-11)
        row = panel_in_page // self.ncol    # Which row (0-3)
        col = panel_in_page % self.ncol     # Which column (0-2)
        
        # Determine if this is the bottom-most panel in its column
        is_bottom_in_column = True
        for check_i in range(i + 1, min(i + (self.nrow - row) * self.ncol, self.nions)):
            check_panel = check_i % self.plotppage
            check_col = check_panel % self.ncol
            if check_col == col:  # Same column, but lower row exists
                is_bottom_in_column = False
                break
        
        # Set x-axis labels and ticks
        if is_bottom_in_column:
            # Show x-tick labels and add x-label
            ax.tick_params(axis='x', labelbottom=True)
            ax.set_xlabel('Velocity [km/s]', fontsize=12)
        else:
            # Hide x-tick labels
            ax.tick_params(axis='x', labelbottom=False)
        
        # Set y-axis label for leftmost column
        if col == 0:  # First column
            ax.set_ylabel('Flux', fontsize=12)
        else:
            ax.set_ylabel('')  # Clear any existing y-label

    
    def plotText(self, flag=1):
        """
        Return text description based on flag value
        
        Parameters:
        -----------
        flag : int, optional
            Flag value (default: 1)
            
        Returns:
        --------
        str
            Text description
        """
        if flag == 1:
            text = 'Detection'       
        elif flag == 0:
            text = 'Non-Detection'      
        elif flag == 2:
            text = 'Blended-detection'
        elif flag == 3:
            text = 'Low-Confidence Detection'
        return text

# Test code for running vStack independently
if __name__ == "__main__":
    from PyQt5.QtWidgets import QApplication
    import sys
    
    print("Testing vStack standalone")
    
    app = QApplication(sys.argv)
    
    # Create a simple dummy parent
    class DummyParent:
        def __init__(self):
            self.canvas = None
            
    parent = DummyParent()
    
    # Create some test data
    wave = np.linspace(1000, 2000, 1000)
    flux = np.random.normal(1.0, 0.1, 1000)
    error = np.random.normal(0.1, 0.01, 1000)
    
    # Create vStack instance
    vs = vStack(parent, wave, flux, error, 'LLS', zabs=0.1)
    
    sys.exit(app.exec_())