# rbcodes/GUIs/multispecviewer/utils.py

def read_line_options():
    """
    Read available line list options from the configuration file.
    
    Returns:
        list: List of line list option names. Returns default options if the
              configuration file cannot be read.
    """
    try:
        # Use pkg_resources to find the configuration file
        from pkg_resources import resource_filename
        config_file = resource_filename('rbcodes.GUIs.multispecviewer', 'line_options.conf')
        
        # Default options in case file is not found
        default_options = ["None", "LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong", "AGN"]
        
        # Try to read the file
        with open(config_file, 'r') as f:
            # Read lines and strip whitespace
            options = [line.strip() for line in f.readlines()]
            
            # Filter out empty lines and comments (lines starting with #)
            options = [opt for opt in options if opt and not opt.startswith('#')]
            
            # If no options were read, use defaults
            if not options:
                print("Warning: No line options found in configuration file, using defaults.")
                return default_options
                
            return options
            
    except Exception as e:
        # If there's any error reading the file, use defaults
        print(f"Error reading line options configuration: {str(e)}")
        print("Using default line options instead.")
        return ["None", "LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong", "AGN"]



def show_help_dialog(parent=None):
    """
    Display a help dialog with keyboard shortcuts and usage information.
    
    Parameters:
    -----------
    parent : QWidget, optional
        Parent widget for the dialog
    """
    from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QTextEdit, 
                                QPushButton, QHBoxLayout, QLabel)
    from PyQt5.QtCore import Qt
    from PyQt5.QtGui import QFont
    import os
    
    # Help text
    help_text = """MULTISPECVIEWER HELP

OVERVIEW:
---------
MultispecViewer is a tool for visualizing and analyzing multiple spectroscopic datasets simultaneously.

INTERFACE:
----------
- Left Widget: Absorber manager for adding, tracking, and plotting absorber systems
- Main Canvas: Main spectral display area with keyboard shortcuts (click here to enable keyboard events)
- Redshift Widget: Shows active redshift for line identification; use Catalog button to add to manager

BASIC NAVIGATION:
----------------
- r: Reset view to original state and clear all lines
- R: Keep spectra but remove all lines/text from canvas
- x: Set left limit (xmin) to cursor position
- X: Set right limit (xmax) to cursor position
- [: Shift view left
- ]: Shift view right
- o: Zoom out (x-axis)
- Y: Manually input y-limits range

DISPLAY CONTROLS:
----------------
- t: Set maximum y-limit to cursor position
- b: Set minimum y-limit to cursor position
- S: Increase smoothing (convolution kernel size)
- U: Decrease smoothing (convolution kernel size)

ANALYSIS TOOLS:
--------------
- v: Open vStack GUI for identifying transitions (press S to exit)
- V: Same as 'v' but allows manual velocity axis selection
- Right-Click: Shows a menu of possible spectral lines at cursor position for identification

QUICK LINE IDENTIFICATION:
-------------------------
- M: Mark MgII doublet (2796, 2803)
- C: Mark CIV doublet (1548, 1550)
- F: Mark FeII multiplet (2600, 2586, 2382)
- 6: Mark OVI doublet (1031, 1037)
- 4: Mark SiIV doublet (1393, 1402)
- 8: Mark NeVIII doublet (778, 770)
- 2: Mark Lyb/Lya
- 1: Mark Lya/Lyb

VSTACK CONTROLS:
---------------
- >: Next page
- <: Previous page
- w: Toggle transition flag between:
   - Detection
   - Blended-detection
   - Low-Confidence Detection
   - Non-Detection
- S: Save transition list and return to main view

FILE MANAGEMENT:
---------------
- Save: Save data in different formats:
   - JSON: Combined data with absorbers, line lists, and metadata
   - CSV: Absorber data only (.csv)
   - TXT: Line list data only (.txt)
- Load: Load previously saved data:
   - Select a JSON file for complete data
   - Select a CSV or TXT file for individual components
   - Choose to append or overwrite existing data

GETTING HELP:
------------
- H: Show this help window
- For complete documentation, see the multispec.md file
"""
    
    # Create dialog
    dialog = QDialog(parent)
    dialog.setWindowTitle("MultispecViewer Help")
    dialog.setMinimumSize(700, 500)
    
    # Create layout
    layout = QVBoxLayout(dialog)
    
    # Create text edit for help content
    text_edit = QTextEdit()
    text_edit.setReadOnly(True)
    text_edit.setPlainText(help_text)
    
    # Set monospace font for better formatting
    font = QFont("Courier New", 10)
    text_edit.setFont(font)
    
    # Add to layout
    layout.addWidget(text_edit)
    
    # Add a label for README link
    readme_label = QLabel("For more detailed documentation:")
    layout.addWidget(readme_label)
    
    # Add buttons for documentation access
    button_layout = QHBoxLayout()

    # Add GitHub button - always available
    github_link = "https://github.com/rongmon/rbcodes/blob/master/docs/GUIs/multispec/multispec.md"
    github_button = QPushButton("View README on GitHub")

    def open_github():
        import webbrowser
        webbrowser.open(github_link)

    github_button.clicked.connect(open_github)
    button_layout.addWidget(github_button)

    # Try to find local README file as well
    try:
        import os
        from pkg_resources import resource_filename
        
        # Check possible local locations
        potential_paths = [
            resource_filename('rbcodes', 'docs/GUIs/multispec/multispec.md'),
            resource_filename('rbcodes.GUIs.multispecviewer', 'multispec.md'),
            resource_filename('rbcodes', 'GUIs/multispecviewer/multispec.md')
        ]
        
        # Find first existing path
        readme_path = None
        for path in potential_paths:
            if os.path.exists(path):
                readme_path = path
                break
        
        # If found locally, add button to open it
        if readme_path:
            local_button = QPushButton("Open Local README")
            
            def open_local_readme():
                import sys
                import subprocess
                
                # Different methods to open file based on OS
                if sys.platform.startswith('darwin'):  # macOS
                    subprocess.call(('open', readme_path))
                elif sys.platform.startswith('win'):   # Windows
                    os.startfile(readme_path)
                else:  # Linux and other Unix-like
                    subprocess.call(('xdg-open', readme_path))
            
            local_button.clicked.connect(open_local_readme)
            button_layout.addWidget(local_button)
            
    except Exception as e:
        print(f"Could not locate local README.md: {str(e)}")

    # Add close button
    close_button = QPushButton("Close")
    close_button.clicked.connect(dialog.accept)
    button_layout.addWidget(close_button)

    # Add button layout to main layout
    layout.addLayout(button_layout)
    
    # Apply dark theme if parent has it
    if parent and hasattr(parent, 'palette'):
        dialog.setPalette(parent.palette())
        
        # Additional styling for dark theme
        dialog.setStyleSheet("""
            QDialog {
                background-color: #353535;
                color: #F2F2F7;
            }
            QTextEdit {
                background-color: #252525;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 6px;
                padding: 8px;
            }
            QLabel {
                color: #F2F2F7;
            }
            QPushButton {
                background-color: #474747;
                color: #F2F2F2;
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050;
            }
            QPushButton:pressed {
                background-color: #2A2A2A;
            }
        """)
    
    # Show dialog
    dialog.exec_()

def reconcile_linelists(input_files, velocity_threshold=20, output_file=None, create_absorber_df=True):
    """
    Reconcile multiple linelist files by merging duplicate entries based on line name and velocity separation.
    
    Parameters:
    -----------
    input_files : list
        List of file paths to linelist files (can be txt, csv, or json format)
    velocity_threshold : float, optional
        Maximum velocity difference in km/s to consider lines as duplicates (default: 20)
    output_file : str, optional
        Path to save the reconciled linelist (if None, doesn't save to file)
    create_absorber_df : bool, optional
        Whether to create an absorber DataFrame from unique redshifts (default: True)
        
    Returns:
    --------
    tuple
        (reconciled_linelist, absorber_df) where reconciled_linelist is a DataFrame with
        merged line entries and absorber_df is a DataFrame of unique absorber systems
    """
    import pandas as pd
    import numpy as np
    import os
    from scipy import constants
    
    # Speed of light in km/s
    c = constants.c / 1000  # Convert from m/s to km/s
    
    # Initialize empty master linelist
    master_linelist = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs'])
    
    # Try to import IO manager for file loading
    try:
        from rbcodes.GUIs.multispecviewer.io_manager import IOManager
        io_manager = IOManager()
        has_io_manager = True
    except ImportError:
        has_io_manager = False
        print("IOManager not available. Using basic file loading methods.")
    
    # Load all input files
    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"Warning: File not found: {file_path}")
            continue
            
        # Determine file format from extension
        ext = os.path.splitext(file_path)[1].lower()
        
        try:
            # Use IOManager if available
            if has_io_manager:
                if ext == '.json':
                    ll_df, _, _, _, _ = io_manager.load_combined_data(file_path)
                else:
                    ll_df, _ = io_manager.load_line_list(file_path)
            # Otherwise use basic loading methods
            else:
                if ext == '.json':
                    # Try to load as JSON
                    import json
                    with open(file_path, 'r') as f:
                        data = json.load(f)
                    
                    if 'line_list' in data:
                        ll_df = pd.DataFrame(data['line_list'])
                    else:
                        # Try different structure
                        ll_df = pd.DataFrame(data)
                elif ext == '.csv':
                    # Load as CSV
                    ll_df = pd.read_csv(file_path)
                else:
                    # Load as text file
                    ll_df = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs'])
                    
                    # Parse the text file
                    with open(file_path, 'r') as f:
                        lines = f.readlines()
                    
                    # Skip header (first two lines)
                    data_lines = lines[2:] if len(lines) > 2 else lines
                    
                    # Parse each line
                    for line in data_lines:
                        if line.strip():  # Skip empty lines
                            try:
                                # Split by whitespace with careful handling
                                parts = line.strip().split()
                                if len(parts) >= 3:
                                    # Extract values
                                    wave_obs = float(parts[-2])
                                    zabs = float(parts[-1])
                                    
                                    # Everything before these two values is the Name
                                    name = ' '.join(parts[:-2])
                                    
                                    # Add to ll_df
                                    new_row = pd.Series({'Name': name, 'Wave_obs': wave_obs, 'Zabs': zabs})
                                    ll_df = pd.concat([ll_df, pd.DataFrame([new_row])], ignore_index=True)
                            except Exception as e:
                                print(f"Error parsing line: {line} - {str(e)}")
                                continue
            
            # Skip if loading failed
            if ll_df is None or ll_df.empty:
                print(f"Warning: No valid data found in {file_path}")
                continue
                
            # Verify required columns
            required_cols = ['Name', 'Wave_obs', 'Zabs']
            if not all(col in ll_df.columns for col in required_cols):
                print(f"Warning: File {file_path} is missing required columns: {required_cols}")
                continue
            
            # Append to master linelist
            master_linelist = pd.concat([master_linelist, ll_df], ignore_index=True)
            print(f"Loaded {len(ll_df)} lines from {file_path}")
            
        except Exception as e:
            print(f"Error loading {file_path}: {str(e)}")
    
    # Check if we have any data to process
    if master_linelist.empty:
        print("No valid linelist data loaded.")
        return master_linelist, None
        
    print(f"Total lines loaded: {len(master_linelist)}")
    
    # Calculate rest wavelength for each line
    master_linelist['Wave_rest'] = master_linelist['Wave_obs'] / (1 + master_linelist['Zabs'])
    
    # Extract base transition names by removing any modifiers
    def extract_base_name(name):
        # Remove annotation markers like [b] or [p]
        base_name = name.split('[')[0].strip()
        return base_name
    
    master_linelist['BaseName'] = master_linelist['Name'].apply(extract_base_name)
    
    # Sort by BaseName and Wave_rest for easier processing
    master_linelist = master_linelist.sort_values(['BaseName', 'Wave_rest']).reset_index(drop=True)
    
    # Group by BaseName
    grouped = master_linelist.groupby('BaseName')
    
    # Initialize reconciled linelist
    reconciled_lines = []
    
    # Process each group
    for name, group in grouped:
        # Skip if only one entry for this name
        if len(group) == 1:
            reconciled_lines.append(group.iloc[0].to_dict())
            continue
            
        # Sort by Wave_rest
        sorted_group = group.sort_values('Wave_rest')
        
        # Initialize cluster with first line
        current_cluster = [sorted_group.iloc[0]]
        
        # Compare each subsequent line
        for i in range(1, len(sorted_group)):
            line = sorted_group.iloc[i]
            last_line = current_cluster[-1]
            
            # Calculate velocity difference in rest frame
            # v = c * (λ2 - λ1) / λ1
            v_diff = c * (line['Wave_rest'] - last_line['Wave_rest']) / last_line['Wave_rest']
            v_diff_abs = abs(v_diff)
            
            # If within threshold, add to current cluster
            if v_diff_abs <= velocity_threshold:
                current_cluster.append(line)
            else:
                # Process current cluster
                if len(current_cluster) > 1:
                    # Create merged entry
                    cluster_df = pd.DataFrame(current_cluster)
                    
                    # Calculate mean redshift
                    mean_z = cluster_df['Zabs'].mean()
                    
                    # Calculate mean rest wavelength
                    mean_rest = cluster_df['Wave_rest'].mean()
                    
                    # Convert back to observed wavelength using mean redshift
                    mean_obs = mean_rest * (1 + mean_z)
                    
                    merged_entry = {
                        'Name': name,
                        'Wave_obs': mean_obs,
                        'Zabs': mean_z,
                        'Wave_rest': mean_rest,
                        'MergedCount': len(current_cluster)
                    }
                    reconciled_lines.append(merged_entry)
                else:
                    # Single entry, no need to merge
                    reconciled_lines.append(current_cluster[0].to_dict())
                
                # Start new cluster
                current_cluster = [line]
        
        # Process final cluster
        if len(current_cluster) > 1:
            # Create merged entry
            cluster_df = pd.DataFrame(current_cluster)
            
            # Calculate mean redshift
            mean_z = cluster_df['Zabs'].mean()
            
            # Calculate mean rest wavelength
            mean_rest = cluster_df['Wave_rest'].mean()
            
            # Convert back to observed wavelength using mean redshift
            mean_obs = mean_rest * (1 + mean_z)
            
            merged_entry = {
                'Name': name,
                'Wave_obs': mean_obs,
                'Zabs': mean_z,
                'Wave_rest': mean_rest,
                'MergedCount': len(current_cluster)
            }
            reconciled_lines.append(merged_entry)
        else:
            # Single entry, no need to merge
            reconciled_lines.append(current_cluster[0].to_dict())
    
    # Create reconciled DataFrame
    reconciled_df = pd.DataFrame(reconciled_lines)
    
    # Clean up DataFrame - remove unnecessary columns
    columns_to_keep = ['Name', 'Wave_obs', 'Zabs']
    if 'MergedCount' in reconciled_df.columns:
        columns_to_keep.append('MergedCount')
    
    reconciled_df = reconciled_df[columns_to_keep]
    
    # Round numerical values
    if 'Wave_obs' in reconciled_df.columns:
        reconciled_df['Wave_obs'] = reconciled_df['Wave_obs'].round(4)
    if 'Zabs' in reconciled_df.columns:
        reconciled_df['Zabs'] = reconciled_df['Zabs'].round(6)
    
    print(f"Reconciled to {len(reconciled_df)} unique lines")
    
    # Create absorber DataFrame if requested
    absorber_df = None
    if create_absorber_df:
        # Calculate redshift tolerance from velocity threshold
        # v/c = Δz/(1+z) -> Δz = (v/c)*(1+z)
        # For small redshifts, we can approximate z_tolerance = v/c
        z_tolerance = velocity_threshold / c
        
        print(f"Using redshift tolerance of {z_tolerance:.6f} (corresponding to {velocity_threshold} km/s)")
        
        # Get unique redshifts from the reconciled linelist
        unique_zabs = []
        
        # Sort by redshift
        sorted_by_z = reconciled_df.sort_values('Zabs').reset_index(drop=True)
        
        # Group similar redshifts
        current_z_group = []  # Will store redshift values
        
        # Add first redshift
        first_z = sorted_by_z['Zabs'].iloc[0]
        current_z_group.append(first_z)
        
        for i in range(1, len(sorted_by_z)):
            current_z = sorted_by_z['Zabs'].iloc[i]
            prev_z = current_z_group[-1]
            
            # Calculate proper velocity difference
            v_diff = c * (current_z - prev_z)/(1 + prev_z)
            v_diff_abs = abs(v_diff)
            
            # If within velocity threshold, add to current group
            if v_diff_abs <= velocity_threshold:
                current_z_group.append(current_z)
            else:
                # Process current group
                if current_z_group:
                    # Calculate mean z
                    mean_z = np.mean(current_z_group)
                    
                    # Add to unique redshifts
                    unique_zabs.append(round(mean_z, 6))
                    
                # Start new group with current z
                current_z_group = [current_z]
        
        # Process final group
        if current_z_group:
            # Calculate mean z
            mean_z = np.mean(current_z_group)
            
            # Add to unique redshifts
            unique_zabs.append(round(mean_z, 6))
        
        # Create initial absorber DataFrame with default values
        absorber_data = []
        for z in unique_zabs:
            absorber_data.append({
                'Zabs': z,
                'LineList': 'LLS',  # Default line list
                'Visible': False    # Default to not visible
            })
        
        absorber_df = pd.DataFrame(absorber_data)
        
        # Now check if any CSV files were provided that might contain absorber information
        csv_files = [f for f in input_files if f.lower().endswith('.csv')]
        
        if csv_files:
            # Try to load absorber information from CSV files
            for csv_file in csv_files:
                try:
                    # Load the CSV
                    csv_data = pd.read_csv(csv_file)
                    
                    # Check if it has required columns for an absorber file
                    if 'Zabs' in csv_data.columns:
                        print(f"Found absorber information in {csv_file}")
                        
                        # Check for LineList column (might be named differently)
                        linelist_col = None
                        for col in csv_data.columns:
                            if col.lower() == 'linelist' or col.lower() == 'list':
                                linelist_col = col
                                break
                        
                        if linelist_col:
                            # Cross-match with absorber_df and update LineList for matching entries
                            for i, row in absorber_df.iterrows():
                                # Find closest matching redshift in CSV data
                                closest_idx = csv_data['Zabs'].sub(row['Zabs']).abs().idxmin()
                                closest_z = csv_data.loc[closest_idx, 'Zabs']
                                
                                # Calculate velocity difference
                                v_diff = c * (closest_z - row['Zabs'])/(1 + row['Zabs'])
                                v_diff_abs = abs(v_diff)
                                
                                # If within velocity threshold, use the LineList from CSV
                                if v_diff_abs <= velocity_threshold:
                                    absorber_df.at[i, 'LineList'] = csv_data.loc[closest_idx, linelist_col]
                                    print(f"Updated LineList for z={row['Zabs']} to {csv_data.loc[closest_idx, linelist_col]}")
                except Exception as e:
                    print(f"Error processing CSV file {csv_file}: {str(e)}")
        
        # Get color list from rb_utility
        try:
            from rbcodes.utils import rb_utility as rt
            clr = rt.rb_set_color()
            colors = list(clr.keys())[1:]  # Skip the first color (usually background)
        except ImportError:
            # Fallback colors if rb_utility not available
            colors = ['white', 'red', 'blue', 'green', 'yellow', 'cyan', 'magenta', 
                    'orange', 'purple', 'pink', 'teal', 'lime', 'brown', 'navy']
        
        # Assign different colors to each absorber, cycling through the color list
        absorber_df['Color'] = [colors[i % len(colors)] for i in range(len(absorber_df))]
        
        print(f"Created {len(absorber_df)} unique absorber systems")
    
    # Save to file if requested
    if output_file:
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        # Determine output format
        output_ext = os.path.splitext(output_file)[1].lower()
        
        if output_ext == '.json':
            # Create combined data structure
            output_data = {
                'line_list': reconciled_df.to_dict(orient='records')
            }
            
            # Add absorber data if available
            if absorber_df is not None and not absorber_df.empty:
                output_data['absorbers'] = absorber_df.to_dict(orient='records')
                
            # Add metadata
            import datetime
            output_data['metadata'] = {
                'creation_date': datetime.datetime.now().isoformat(),
                'reconciliation_info': {
                    'velocity_threshold': velocity_threshold,
                    'input_files': [os.path.basename(f) for f in input_files],
                    'original_line_count': len(master_linelist),
                    'reconciled_line_count': len(reconciled_df)
                }
            }
            
            # Save to JSON
            import json
            with open(output_file, 'w') as f:
                json.dump(output_data, f, indent=2)
                
            print(f"Saved reconciled data to {output_file}")
        else:
            # Warn about extension but save as JSON anyway
            print(f"Warning: Recommended file extension is .json. Saving as JSON format.")
            
            output_file_json = os.path.splitext(output_file)[0] + '.json'
            
            # Create combined data structure
            output_data = {
                'line_list': reconciled_df.to_dict(orient='records')
            }
            
            # Add absorber data if available
            if absorber_df is not None and not absorber_df.empty:
                output_data['absorbers'] = absorber_df.to_dict(orient='records')
                
            # Add metadata
            import datetime
            output_data['metadata'] = {
                'creation_date': datetime.datetime.now().isoformat(),
                'reconciliation_info': {
                    'velocity_threshold': velocity_threshold,
                    'input_files': [os.path.basename(f) for f in input_files],
                    'original_line_count': len(master_linelist),
                    'reconciled_line_count': len(reconciled_df)
                }
            }
            
            # Save to JSON
            import json
            with open(output_file_json, 'w') as f:
                json.dump(output_data, f, indent=2)
                
            print(f"Saved reconciled data to {output_file_json}")
    
    return reconciled_df, absorber_df