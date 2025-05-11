# io_manager.py

import os
import json
import datetime
import pandas as pd
import platform
import getpass
import traceback
from pathlib import Path

class IOManager:
    """
    Singleton class to handle all file I/O operations for MultispecViewer.
    Provides methods for loading and saving various file formats including
    JSON for combined data, and traditional text/CSV formats for backward compatibility.
    """
    _instance = None
    
    def __new__(cls):
        """Ensure only one instance of IOManager exists (singleton pattern)"""
        if cls._instance is None:
            cls._instance = super(IOManager, cls).__new__(cls)
            cls._instance._initialize()
        return cls._instance
    
    def _initialize(self):
        """Initialize the IOManager with default settings"""
        # Store last used directories for different file types
        self.last_directories = {
            'fits': os.getcwd(),
            'line_list': os.getcwd(),
            'absorbers': os.getcwd(),
            'json': os.getcwd(),
        }
        
        # Version info for file formats
        self.current_version = "1.0.0"
        
        # Reference to message box (to be set by MainWindow)
        self.message_box = None
    
    def set_message_box(self, message_box):
        """Set the message box reference for displaying messages"""
        self.message_box = message_box
    
    def show_message(self, message, color="#FFFFFF"):
        """Display a message in the message box if available, otherwise print"""
        if self.message_box:
            self.message_box.on_sent_message(message, color)
        else:
            print(message)
    
    # ===== FITS File Operations =====
    
    def load_fits_files(self, file_paths=None):
        """
        Load FITS files into XSpectrum1D objects.
        
        Args:
            file_paths (list, optional): List of file paths to load. If None, open file dialog.
            
        Returns:
            list: List of XSpectrum1D objects
        """
        try:
            # Import here to avoid circular imports
            from linetools.spectra.xspectrum1d import XSpectrum1D
            import numpy as np
            
            if file_paths is None:
                # This would be handled by the calling code that has access to the UI
                return [], "No files specified"
            
            spectra = []
            loaded_files = []
            errors = []
            
            for file_path in file_paths:
                try:
                    if not os.path.exists(file_path):
                        errors.append(f"File not found: {file_path}")
                        continue
                        
                    # Update last directory
                    self.last_directories['fits'] = os.path.dirname(file_path)
                    
                    # Load the spectrum
                    temp_spec = XSpectrum1D.from_file(file_path)
                    wave = temp_spec.wavelength.value
                    flux = temp_spec.flux.value
                    
                    # Check if the original spectrum has an error array
                    if temp_spec.sig_is_set:
                        # Use the existing error array
                        sig = temp_spec.sig.value
                        spec = XSpectrum1D.from_tuple((wave, flux, sig), verbose=False)
                    else:
                        # Create error array as 5% of flux values
                        sig_array = 0.05 * flux
                        spec = XSpectrum1D.from_tuple((wave, flux, sig_array), verbose=False)
                        self.show_message(
                            f"No error spectrum found for {os.path.basename(file_path)}. " 
                            f"Using 5% of flux values as error.", "#FFA500")
                    
                    # Set the filename attribute to match the original
                    spec.filename = file_path
                    
                    spectra.append(spec)
                    loaded_files.append(os.path.basename(file_path))
                    
                except Exception as e:
                    errors.append(f"Error loading {os.path.basename(file_path)}: {str(e)}")
            
            # Show summary message
            if spectra:
                self.show_message(f"Successfully loaded {len(spectra)} files: {', '.join(loaded_files)}", "#008000")
            
            if errors:
                error_msg = "\n".join(errors)
                self.show_message(f"Errors encountered:\n{error_msg}", "#FF0000")
            
            return spectra, None if not errors else "\n".join(errors)
            
        except Exception as e:
            error_message = f"Error in load_fits_files: {str(e)}"
            self.show_message(error_message, "#FF0000")
            return [], error_message
    
    # ===== Line List Operations =====   
    def save_line_list(self, line_list, file_path=None, format=None):
        """
        Save a line list DataFrame to a file.
        
        Args:
            line_list (pd.DataFrame): DataFrame containing the line list
            file_path (str, optional): Path to save to. If None, use the last directory.
            format (str, optional): Format to save in ('txt', 'csv', 'json'). 
                                    If None, determine from file extension.
        
        Returns:
            bool: Success status
            str: Error message if any
        """
        try:
            if line_list is None or line_list.empty:
                return False, "No line list data to save"
            
            if file_path is None:
                return False, "No file path specified"
            
            # Determine format if not specified
            if format is None:
                ext = os.path.splitext(file_path)[1].lower()
                if ext == '.txt':
                    format = 'txt'
                elif ext == '.csv':
                    format = 'csv'
                elif ext == '.json':
                    format = 'json'
                else:
                    format = 'txt'  # Default to txt
            
            # Update last directory
            self.last_directories['line_list'] = os.path.dirname(file_path)
            
            # Save based on format
            if format == 'txt':
                with open(file_path, 'w') as f:
                    # Write header
                    f.write(f"{'Name':<30} {'Wave_obs':<15} {'Zabs':<10}\n")
                    f.write("-" * 55 + "\n")
                    
                    # Write each line
                    for _, row in line_list.iterrows():
                        name = str(row['Name'])
                        wave_obs = f"{row['Wave_obs']:.4f}"
                        zabs = f"{row['Zabs']:.6f}"
                        f.write(f"{name:<30} {wave_obs:<15} {zabs:<10}\n")
                
                self.show_message(f"Saved line list to: {file_path}", "#008000")
                return True, None
                
            elif format == 'csv':
                line_list.to_csv(file_path, index=False)
                self.show_message(f"Saved line list to: {file_path}", "#008000")
                return True, None
                
            elif format == 'json':
                # Convert to dict for JSON serialization
                line_list_dict = line_list.to_dict(orient='records')
                with open(file_path, 'w') as f:
                    json.dump({'line_list': line_list_dict}, f, indent=2)
                
                self.show_message(f"Saved line list to: {file_path}", "#008000")
                return True, None
            
            else:
                return False, f"Unknown format: {format}"
                
        except Exception as e:
            error_message = f"Error saving line list: {str(e)}"
            self.show_message(error_message, "#FF0000")
            return False, error_message
    
    def load_line_list(self, file_path=None):
        """
        Load a line list from a file.
        
        Args:
            file_path (str, optional): Path to the file to load
            
        Returns:
            pd.DataFrame: The loaded line list
            str: Error message if any
        """
        try:
            if file_path is None:
                return None, "No file path specified"
            
            if not os.path.exists(file_path):
                return None, f"File not found: {file_path}"
            
            # Update last directory
            self.last_directories['line_list'] = os.path.dirname(file_path)
            
            # Determine format from extension
            ext = os.path.splitext(file_path)[1].lower()
            
            # Load based on format
            if ext == '.json':
                with open(file_path, 'r') as f:
                    data = json.load(f)
                
                if 'line_list' in data:
                    line_list = pd.DataFrame(data['line_list'])
                    self.show_message(f"Loaded {len(line_list)} line identifications from JSON", "#008000")
                    return line_list, None
                else:
                    # Try to load as a simple JSON array
                    try:
                        line_list = pd.DataFrame(data)
                        self.show_message(f"Loaded {len(line_list)} line identifications from JSON", "#008000")
                        return line_list, None
                    except:
                        return None, "Invalid JSON format: missing 'line_list' key"
            
            elif ext == '.csv':
                try:
                    line_list = pd.read_csv(file_path)
                    # Verify required columns
                    required_cols = ['Name', 'Wave_obs', 'Zabs']
                    if not all(col in line_list.columns for col in required_cols):
                        missing = [col for col in required_cols if col not in line_list.columns]
                        return None, f"CSV file missing required columns: {', '.join(missing)}"
                    
                    self.show_message(f"Loaded {len(line_list)} line identifications from CSV", "#008000")
                    return line_list, None
                except Exception as e:
                    return None, f"Error loading CSV: {str(e)}"
            
            elif ext == '.txt':
                try:
                    # Initialize line_list DataFrame
                    line_list = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs'])
                    
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
                                    wave_index = -2  # Assume Wave_obs is second to last
                                    zabs_index = -1  # Assume Zabs is last
                                    
                                    # Extract the last two values as Wave_obs and Zabs
                                    wave_obs = float(parts[wave_index])
                                    zabs = float(parts[zabs_index])
                                    
                                    # Everything before these two values is the Name
                                    name = ' '.join(parts[:wave_index])
                                    
                                    # Add to line_list
                                    new_row = pd.Series({'Name': name, 'Wave_obs': wave_obs, 'Zabs': zabs})
                                    line_list = pd.concat([line_list, pd.DataFrame([new_row])], ignore_index=True)
                            except ValueError:
                                # Try another approach
                                try:
                                    import re
                                    columns = re.split(r'\s{2,}', line.strip())
                                    if len(columns) >= 3:
                                        name = columns[0].strip()
                                        wave_obs = float(columns[1].strip())
                                        zabs = float(columns[2].strip())
                                        
                                        # Add to line_list
                                        new_row = pd.Series({'Name': name, 'Wave_obs': wave_obs, 'Zabs': zabs})
                                        line_list = pd.concat([line_list, pd.DataFrame([new_row])], ignore_index=True)
                                except Exception as e:
                                    print(f"Could not parse line: {line} - {str(e)}")
                    
                    self.show_message(f"Loaded {len(line_list)} line identifications from text file", "#008000")
                    return line_list, None
                except Exception as e:
                    return None, f"Error parsing text file: {str(e)}"
            
            else:
                return None, f"Unsupported file extension: {ext}"
        
        except Exception as e:
            error_message = f"Error loading line list: {str(e)}"
            self.show_message(error_message, "#FF0000")
            return None, error_message
    
    # ===== Absorber Operations ===== 
    def save_absorbers(self, absorbers_df, file_path=None, format=None):
        """
        Save absorbers DataFrame to a file.
        
        Args:
            absorbers_df (pd.DataFrame): DataFrame containing absorber data
            file_path (str, optional): Path to save to. If None, use the last directory.
            format (str, optional): Format to save in ('csv', 'json'). 
                                     If None, determine from file extension.
        
        Returns:
            bool: Success status
            str: Error message if any
        """
        try:
            if absorbers_df is None or absorbers_df.empty:
                return False, "No absorber data to save"
            
            if file_path is None:
                return False, "No file path specified"
            
            # Determine format if not specified
            if format is None:
                ext = os.path.splitext(file_path)[1].lower()
                if ext == '.csv':
                    format = 'csv'
                elif ext == '.json':
                    format = 'json'
                else:
                    format = 'csv'  # Default to csv
            
            # Update last directory
            self.last_directories['absorbers'] = os.path.dirname(file_path)
            
            # Save based on format
            if format == 'csv':
                # Select only the columns we want to save
                required_cols = ['Zabs', 'LineList', 'Color']
                if all(col in absorbers_df.columns for col in required_cols):
                    save_df = absorbers_df[required_cols]
                else:
                    save_df = absorbers_df
                
                save_df.to_csv(file_path, index=False)
                self.show_message(f"Saved absorber data to: {file_path}", "#008000")
                return True, None
                
            elif format == 'json':
                # Convert to dict for JSON serialization
                absorbers_dict = absorbers_df.to_dict(orient='records')
                with open(file_path, 'w') as f:
                    json.dump({'absorbers': absorbers_dict}, f, indent=2)
                
                self.show_message(f"Saved absorber data to: {file_path}", "#008000")
                return True, None
            
            else:
                return False, f"Unknown format: {format}"
                
        except Exception as e:
            error_message = f"Error saving absorber data: {str(e)}"
            self.show_message(error_message, "#FF0000")
            return False, error_message
    
    def load_absorbers(self, file_path=None):
        """
        Load absorber data from a file.
        
        Args:
            file_path (str, optional): Path to the file to load
            
        Returns:
            pd.DataFrame: The loaded absorber data
            str: Error message if any
        """
        try:
            if file_path is None:
                return None, "No file path specified"
            
            if not os.path.exists(file_path):
                return None, f"File not found: {file_path}"
            
            # Update last directory
            self.last_directories['absorbers'] = os.path.dirname(file_path)
            
            # Determine format from extension
            ext = os.path.splitext(file_path)[1].lower()
            
            # Load based on format
            if ext == '.json':
                with open(file_path, 'r') as f:
                    data = json.load(f)
                
                if 'absorbers' in data:
                    absorbers_df = pd.DataFrame(data['absorbers'])
                    self.show_message(f"Loaded {len(absorbers_df)} absorber systems from JSON", "#008000")
                    return absorbers_df, None
                else:
                    # Try to load as a simple JSON array
                    try:
                        absorbers_df = pd.DataFrame(data)
                        self.show_message(f"Loaded {len(absorbers_df)} absorber systems from JSON", "#008000")
                        return absorbers_df, None
                    except:
                        return None, "Invalid JSON format: missing 'absorbers' key"
            
            elif ext == '.csv':
                try:
                    absorbers_df = pd.read_csv(file_path)
                    # Verify required columns
                    required_cols = ['Zabs', 'LineList', 'Color']
                    if not all(col in absorbers_df.columns for col in required_cols):
                        missing = [col for col in required_cols if col not in absorbers_df.columns]
                        return None, f"CSV file missing required columns: {', '.join(missing)}"
                    
                    self.show_message(f"Loaded {len(absorbers_df)} absorber systems from CSV", "#008000")
                    return absorbers_df, None
                except Exception as e:
                    return None, f"Error loading CSV: {str(e)}"
            
            else:
                return None, f"Unsupported file extension: {ext}"
        
        except Exception as e:
            error_message = f"Error loading absorber data: {str(e)}"
            self.show_message(error_message, "#FF0000")
            return None, error_message
    
    # ===== Combined JSON Operations =====
    
    def save_combined_data(self, line_list, absorbers_df, spectrum_files=None, file_path=None, user_comment=None, metadata=None):
        """
        Save line list, absorber data, and metadata to a single JSON file.
        
        Args:
            line_list (pd.DataFrame): Line list DataFrame
            absorbers_df (pd.DataFrame): Absorbers DataFrame
            spectrum_files (list, optional): List of spectrum filenames
            file_path (str, optional): Path to save to
            user_comment (str, optional): User-provided comment about the data
            metadata (dict, optional): Additional metadata to include
            
        Returns:
            bool: Success status
            str: Error message if any
        """
        try:
            if file_path is None:
                return False, "No file path specified"
                
            # Add .json extension if not present
            if not file_path.lower().endswith('.json'):
                file_path += '.json'
            
            # Update last directory
            self.last_directories['json'] = os.path.dirname(file_path)
            
            # Build the combined data structure
            combined_data = {
                # Line list data
                'line_list': line_list.to_dict(orient='records') if line_list is not None and not line_list.empty else [],
                
                # Absorber data
                'absorbers': absorbers_df.to_dict(orient='records') if absorbers_df is not None and not absorbers_df.empty else [],
                
                # Spectrum files
                'spectrum_files': spectrum_files or [],
                
                # Metadata
                'metadata': {
                    'creation_date': datetime.datetime.now().isoformat(),
                    'version': self.current_version,
                    'application_name': 'MultispecViewer',
                    'user_comment': user_comment or '',
                    'system_info': {
                        'platform': platform.platform(),
                        'python_version': platform.python_version(),
                        'username': getpass.getuser(),
                    }
                }
            }
            
            # Add any additional metadata
            if metadata:
                combined_data['metadata'].update(metadata)
            
            # Save to file
            with open(file_path, 'w') as f:
                json.dump(combined_data, f, indent=2)
            
            self.show_message(f"Saved combined data to: {file_path}", "#008000")
            return True, None
            
        except Exception as e:
            error_message = f"Error saving combined data: {str(e)}"
            self.show_message(error_message, "#FF0000")
            traceback.print_exc()
            return False, error_message
    
    def load_combined_data(self, file_path=None):
        """
        Load combined data from a JSON file.
        
        Args:
            file_path (str, optional): Path to the file to load
            
        Returns:
            tuple: (line_list, absorbers_df, spectrum_files, metadata, error_message)
                   Any of these may be None if not found or if an error occurs
        """
        try:
            if file_path is None:
                return None, None, None, None, "No file path specified"
            
            if not os.path.exists(file_path):
                return None, None, None, None, f"File not found: {file_path}"
            
            # Update last directory
            self.last_directories['json'] = os.path.dirname(file_path)
            
            # Load the file
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Extract line list
            line_list = None
            if 'line_list' in data and data['line_list']:
                line_list = pd.DataFrame(data['line_list'])
            
            # Extract absorbers
            absorbers_df = None
            if 'absorbers' in data and data['absorbers']:
                absorbers_df = pd.DataFrame(data['absorbers'])
            
            # Extract spectrum files
            spectrum_files = data.get('spectrum_files', [])
            
            # Extract metadata
            metadata = data.get('metadata', {})
            
            # Check for required data
            if line_list is None and absorbers_df is None:
                return None, None, spectrum_files, metadata, "No line list or absorber data found in file"
            
            # Show success message
            message = "Loaded "
            if line_list is not None:
                message += f"{len(line_list)} line identifications"
            if absorbers_df is not None:
                if line_list is not None:
                    message += " and "
                message += f"{len(absorbers_df)} absorber systems"
            self.show_message(message, "#008000")
            
            # Display user comment if available
            user_comment = metadata.get('user_comment', '')
            if user_comment:
                self.show_message(f"Comment: {user_comment}", "#FFA500")
            
            return line_list, absorbers_df, spectrum_files, metadata, None
            
        except Exception as e:
            error_message = f"Error loading combined data: {str(e)}"
            self.show_message(error_message, "#FF0000")
            traceback.print_exc()
            return None, None, None, None, error_message
    
    # ===== Conversion Utilities =====
    
    def convert_text_to_json(self, line_list_path=None, absorbers_path=None, output_path=None, 
                             user_comment=None, spectrum_files=None, metadata=None):
        """
        Convert text/CSV files to a combined JSON file.
        
        Args:
            line_list_path (str, optional): Path to line list text/CSV file
            absorbers_path (str, optional): Path to absorbers CSV file
            output_path (str, optional): Path for the output JSON file
            user_comment (str, optional): User-provided comment
            spectrum_files (list, optional): List of spectrum filenames
            metadata (dict, optional): Additional metadata
            
        Returns:
            bool: Success status
            str: Error message if any
        """
        try:
            # Load line list if provided
            line_list = None
            if line_list_path:
                line_list, error = self.load_line_list(line_list_path)
                if error:
                    return False, f"Error loading line list: {error}"
            
            # Load absorbers if provided
            absorbers_df = None
            if absorbers_path:
                absorbers_df, error = self.load_absorbers(absorbers_path)
                if error:
                    return False, f"Error loading absorbers: {error}"
            
            # Verify we have at least one data source
            if line_list is None and absorbers_df is None:
                return False, "No data provided for conversion"
            
            # Generate output path if not provided
            if output_path is None:
                if line_list_path:
                    base_path = os.path.splitext(line_list_path)[0]
                else:
                    base_path = os.path.splitext(absorbers_path)[0]
                output_path = f"{base_path}_combined.json"
            
            # Add additional metadata about the source files
            if metadata is None:
                metadata = {}
            
            metadata['conversion_info'] = {
                'line_list_source': os.path.basename(line_list_path) if line_list_path else None,
                'absorbers_source': os.path.basename(absorbers_path) if absorbers_path else None,
                'conversion_date': datetime.datetime.now().isoformat()
            }
            
            # Save combined data
            result, error = self.save_combined_data(
                line_list, absorbers_df, spectrum_files, output_path, user_comment, metadata
            )
            
            if result:
                self.show_message(f"Successfully converted to JSON: {output_path}", "#008000")
                return True, None
            else:
                return False, f"Error saving combined data: {error}"
                
        except Exception as e:
            error_message = f"Error converting to JSON: {str(e)}"
            self.show_message(error_message, "#FF0000")
            return False, error_message
    
    def convert_json_to_text(self, json_path=None, output_dir=None):
        """
        Convert a combined JSON file to separate text/CSV files.
        
        Args:
            json_path (str, optional): Path to the JSON file
            output_dir (str, optional): Directory to save the output files
            
        Returns:
            tuple: (line_list_path, absorbers_path, metadata_path, error_message)
                   Any of these may be None if not applicable or if an error occurs
        """
        try:
            if not json_path:
                return None, None, None, "No JSON file specified"
            
            # Load the JSON file
            line_list, absorbers_df, spectrum_files, metadata, error = self.load_combined_data(json_path)
            
            if error:
                return None, None, None, f"Error loading JSON: {error}"
            
            # Set output directory
            if output_dir is None:
                output_dir = os.path.dirname(json_path)
            
            # Create output directory if it doesn't exist
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate base name for output files
            base_name = os.path.splitext(os.path.basename(json_path))[0]
            
            # Save line list if available
            line_list_path = None
            if line_list is not None and not line_list.empty:
                line_list_path = os.path.join(output_dir, f"{base_name}_lines.txt")
                success, error = self.save_line_list(line_list, line_list_path, 'txt')
                if not success:
                    return None, None, None, f"Error saving line list: {error}"
            
            # Save absorbers if available
            absorbers_path = None
            if absorbers_df is not None and not absorbers_df.empty:
                absorbers_path = os.path.join(output_dir, f"{base_name}_absorbers.csv")
                success, error = self.save_absorbers(absorbers_df, absorbers_path, 'csv')
                if not success:
                    return line_list_path, None, None, f"Error saving absorbers: {error}"
            
            # Save metadata to a separate file
            metadata_path = os.path.join(output_dir, f"{base_name}_metadata.txt")
            try:
                with open(metadata_path, 'w') as f:
                    f.write(f"MultispecViewer Data Export\n")
                    f.write(f"=========================\n\n")
                    
                    # Write file information
                    f.write(f"Exported from: {json_path}\n")
                    f.write(f"Export date: {datetime.datetime.now().isoformat()}\n\n")
                    
                    # Write metadata
                    f.write(f"Original Metadata\n")
                    f.write(f"----------------\n")
                    for key, value in metadata.items():
                        if isinstance(value, dict):
                            f.write(f"{key}:\n")
                            for subkey, subvalue in value.items():
                                f.write(f"  {subkey}: {subvalue}\n")
                        else:
                            f.write(f"{key}: {value}\n")
                    f.write("\n")
                    
                    # Write spectrum files if available
                    if spectrum_files:
                        f.write(f"Spectrum Files\n")
                        f.write(f"-------------\n")
                        for i, filename in enumerate(spectrum_files):
                            f.write(f"{i+1}. {filename}\n")
                        f.write("\n")
                    
                    # Write data summary
                    f.write(f"Data Summary\n")
                    f.write(f"-----------\n")
                    if line_list is not None:
                        f.write(f"Line identifications: {len(line_list)}\n")
                        f.write(f"Line list saved to: {line_list_path}\n")
                    if absorbers_df is not None:
                        f.write(f"Absorber systems: {len(absorbers_df)}\n")
                        f.write(f"Absorbers saved to: {absorbers_path}\n")
            except Exception as e:
                self.show_message(f"Error saving metadata: {str(e)}", "#FFA500")
                metadata_path = None
            
            # Show success message
            self.show_message(f"Successfully exported to text/CSV files", "#008000")
            
            return line_list_path, absorbers_path, metadata_path, None
            
        except Exception as e:
            error_message = f"Error converting to text: {str(e)}"
            self.show_message(error_message, "#FF0000")
            return None, None, None, error_message
    
    # ===== User Comment Dialog =====
    
    def get_user_comment_dialog(self, parent=None):
            """
            Create a dialog for entering metadata when saving files.
            
            Args:
                parent (QWidget, optional): Parent widget for the dialog
                
            Returns:
                tuple: (comment, metadata_dict, accepted)
                       comment: string with user's comment
                       metadata_dict: dictionary with additional metadata
                       accepted: boolean indicating if dialog was accepted
            """
            try:
                from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QLabel, QTextEdit, 
                                             QLineEdit, QDialogButtonBox, QFormLayout,
                                             QWidget)
                from PyQt5.QtCore import Qt
                
                dialog = QDialog(parent)
                dialog.setWindowTitle("Add Metadata")
                dialog.setMinimumWidth(400)
                dialog.setMinimumHeight(300)
                
                # Apply Fusion style for consistent appearance
                dialog.setStyle(parent.style() if parent else None)
                
                # Create layout
                layout = QVBoxLayout(dialog)
                
                # Add form for metadata fields
                form_layout = QFormLayout()
                
                # Add title field
                title_edit = QLineEdit()
                form_layout.addRow("Title:", title_edit)
                
                # Add author field
                author_edit = QLineEdit()
                form_layout.addRow("Author:", author_edit)
                
                # Add target field
                target_edit = QLineEdit()
                form_layout.addRow("Target Object:", target_edit)
                
                layout.addLayout(form_layout)
                
                # Add comment label
                comment_label = QLabel("Comments (optional):")
                layout.addWidget(comment_label)
                
                # Add comment text area
                comment_edit = QTextEdit()
                comment_edit.setPlaceholderText("Enter any additional notes about this data...")
                layout.addWidget(comment_edit)
                
                # Add button box
                button_box = QDialogButtonBox(
                    QDialogButtonBox.Save | QDialogButtonBox.Cancel)
                button_box.accepted.connect(dialog.accept)
                button_box.rejected.connect(dialog.reject)
                layout.addWidget(button_box)
                
                # Style the dialog for dark theme
                if parent:
                    dialog.setStyleSheet("""
                        QDialog {
                            background-color: #303030;
                            color: #F0F0F0;
                        }
                        QLabel {
                            color: #F0F0F0;
                        }
                        QTextEdit {
                            background-color: #404040;
                            color: #F0F0F0;
                            border: 1px solid #505050;
                            border-radius: 4px;
                        }
                        QLineEdit {
                            background-color: #404040;
                            color: #F0F0F0;
                            border: 1px solid #505050;
                            border-radius: 4px;
                            padding: 2px 4px;
                        }
                        QPushButton {
                            background-color: #404040;
                            color: #F0F0F0;
                            border: 1px solid #505050;
                            border-radius: 4px;
                            padding: 4px 8px;
                        }
                        QPushButton:hover {
                            background-color: #505050;
                        }
                        QPushButton:pressed {
                            background-color: #606060;
                        }
                    """)
                
                # Show the dialog
                result = dialog.exec_()
                
                if result == QDialog.Accepted:
                    # Get values
                    comment = comment_edit.toPlainText()
                    
                    # Create metadata dictionary
                    metadata = {
                        'title': title_edit.text(),
                        'author': author_edit.text(),
                        'target': target_edit.text()
                    }
                    
                    # Remove empty entries
                    metadata = {k: v for k, v in metadata.items() if v}
                    
                    return comment, metadata, True
                else:
                    return "", {}, False
                    
            except Exception as e:
                self.show_message(f"Error creating comment dialog: {str(e)}", "#FF0000")
                return "", {}, False
        
    # ===== Integration Methods for MainWindow =====
    
    def integrated_save_data(self, parent, line_list=None, absorbers_df=None, spectrum_filenames=None):
        """
        Integrated method for saving data with format selection.
        Shows file dialog and optional metadata dialog.
        
        Args:
            parent (QWidget): Parent widget for dialogs
            line_list (pd.DataFrame, optional): Line list data
            absorbers_df (pd.DataFrame, optional): Absorber data
            spectrum_filenames (list, optional): List of spectrum filenames
            
        Returns:
            bool: Success status
        """
        try:
            from PyQt5.QtWidgets import QFileDialog
            
            # Check if we have data to save
            if (line_list is None or line_list.empty) and (absorbers_df is None or absorbers_df.empty):
                self.show_message("No data to save", "#FFA500")
                return False
            
            # Create file dialog
            options = QFileDialog.Options()
            if line_list is not None and not line_list.empty and absorbers_df is not None and not absorbers_df.empty:
                # Both data types available
                file_filter = "JSON Files (*.json);;Line List (*.txt);;Absorbers (*.csv)"
                suggestion = "data.json"
            elif line_list is not None and not line_list.empty:
                # Only line list
                file_filter = "Line List (*.txt);;JSON Files (*.json)"
                suggestion = "line_list.txt"
            else:
                # Only absorbers
                file_filter = "Absorbers (*.csv);;JSON Files (*.json)"
                suggestion = "absorbers.csv"
            
            # Get initial directory
            initial_dir = self.last_directories.get('fits', os.getcwd())
            if not os.path.exists(initial_dir):
                initial_dir = os.getcwd()
            
            # Show save dialog
            file_path, selected_filter = QFileDialog.getSaveFileName(
                parent, "Save Data", os.path.join(initial_dir, suggestion),
                file_filter, options=options
            )
            
            if not file_path:
                return False  # User cancelled
            
            # Determine format from selected filter
            if "JSON" in selected_filter:
                format_type = "json"
                # Show metadata dialog for JSON
                comment, additional_metadata, accepted = self.get_user_comment_dialog(parent)
                if not accepted:
                    return False  # User cancelled
                
                # Ensure .json extension
                if not file_path.lower().endswith('.json'):
                    file_path += '.json'
                
                # Save combined data
                success, error = self.save_combined_data(
                    line_list, absorbers_df, spectrum_filenames, file_path, comment, additional_metadata
                )
                
                if not success:
                    self.show_message(f"Error saving data: {error}", "#FF0000")
                    return False
                
                return True
                
            elif "Line List" in selected_filter:
                # Ensure .txt extension
                if not file_path.lower().endswith('.txt'):
                    file_path += '.txt'
                
                # Save line list
                success, error = self.save_line_list(line_list, file_path, 'txt')
                
                if not success:
                    self.show_message(f"Error saving line list: {error}", "#FF0000")
                    return False
                
                # If we also have absorbers, offer to save them
                if absorbers_df is not None and not absorbers_df.empty:
                    absorbers_path = os.path.splitext(file_path)[0] + "_absorbers.csv"
                    success, error = self.save_absorbers(absorbers_df, absorbers_path, 'csv')
                    
                    if not success:
                        self.show_message(f"Error saving absorbers: {error}", "#FF0000")
                
                return True
                
            elif "Absorbers" in selected_filter:
                # Ensure .csv extension
                if not file_path.lower().endswith('.csv'):
                    file_path += '.csv'
                
                # Save absorbers
                success, error = self.save_absorbers(absorbers_df, file_path, 'csv')
                
                if not success:
                    self.show_message(f"Error saving absorbers: {error}", "#FF0000")
                    return False
                
                # If we also have line list, offer to save it
                if line_list is not None and not line_list.empty:
                    line_list_path = os.path.splitext(file_path)[0] + "_lines.txt"
                    success, error = self.save_line_list(line_list, line_list_path, 'txt')
                    
                    if not success:
                        self.show_message(f"Error saving line list: {error}", "#FF0000")
                
                return True
            
            else:
                self.show_message("Unknown file format selected", "#FF0000")
                return False
                
        except Exception as e:
            self.show_message(f"Error saving data: {str(e)}", "#FF0000")
            traceback.print_exc()
            return False
    
    def integrated_load_data(self, parent):
        """
        Integrated method for loading data with format detection.
        Shows file dialog and loads data accordingly.
        
        Args:
            parent (QWidget): Parent widget for dialogs
            
        Returns:
            tuple: (line_list, absorbers_df, error_message)
        """
        try:
            from PyQt5.QtWidgets import QFileDialog
            
            # Create file dialog
            options = QFileDialog.Options()
            file_filter = "All Supported Files (*.json *.txt *.csv);;JSON Files (*.json);;Line List (*.txt);;Absorbers (*.csv);;All Files (*)"
            
            # Get initial directory
            initial_dir = self.last_directories.get('fits', os.getcwd())
            if not os.path.exists(initial_dir):
                # Fall back to other directories in order of preference
                for dir_type in ['json', 'line_list', 'absorbers']:
                    if dir_type in self.last_directories and os.path.exists(self.last_directories[dir_type]):
                        initial_dir = self.last_directories[dir_type]
                        break
                else:
                    initial_dir = os.getcwd()
            
            # Show open dialog
            file_path, selected_filter = QFileDialog.getOpenFileName(
                parent, "Load Data", initial_dir,
                file_filter, options=options
            )
            
            if not file_path:
                return None, None, "Operation cancelled"  # User cancelled
            
            # Determine format from file extension
            ext = os.path.splitext(file_path)[1].lower()
            
            if ext == '.json':
                # Load JSON data
                line_list, absorbers_df, spectrum_files, metadata, error = self.load_combined_data(file_path)
                
                if error:
                    self.show_message(f"Error loading data: {error}", "#FF0000")
                    return None, None, error
                
                # Display metadata summary if available
                if metadata:
                    # Show title if available
                    if 'title' in metadata and metadata['title']:
                        self.show_message(f"Title: {metadata['title']}", "#8AB4F8")
                    
                    # Show creation date if available
                    if 'creation_date' in metadata:
                        try:
                            date = datetime.datetime.fromisoformat(metadata['creation_date'])
                            self.show_message(f"Created: {date.strftime('%Y-%m-%d %H:%M')}", "#8AB4F8")
                        except:
                            pass
                
                return line_list, absorbers_df, None
                
            elif ext == '.txt':
                # Load line list
                line_list, error = self.load_line_list(file_path)
                
                if error:
                    self.show_message(f"Error loading line list: {error}", "#FF0000")
                    return None, None, error
                
                # Check for matching absorbers file
                absorbers_path = os.path.splitext(file_path)[0] + "_absorbers.csv"
                if os.path.exists(absorbers_path):
                    absorbers_df, error = self.load_absorbers(absorbers_path)
                    
                    if error:
                        self.show_message(f"Error loading absorbers: {error}", "#FFA500")
                        absorbers_df = None
                else:
                    absorbers_df = None
                
                return line_list, absorbers_df, None
                
            elif ext == '.csv':
                # Try loading as absorbers first
                absorbers_df, error = self.load_absorbers(file_path)
                
                if not error:
                    # Check for matching line list file
                    line_list_path = os.path.splitext(file_path)[0] + "_lines.txt"
                    if os.path.exists(line_list_path):
                        line_list, error = self.load_line_list(line_list_path)
                        
                        if error:
                            self.show_message(f"Error loading line list: {error}", "#FFA500")
                            line_list = None
                    else:
                        line_list = None
                    
                    return line_list, absorbers_df, None
                
                # If failed as absorbers, try as line list
                line_list, error2 = self.load_line_list(file_path)
                
                if not error2:
                    return line_list, None, None
                
                # If both failed, return the first error
                self.show_message(f"Error loading file: {error}", "#FF0000")
                return None, None, error
            
            else:
                error_msg = f"Unsupported file extension: {ext}"
                self.show_message(error_msg, "#FF0000")
                return None, None, error_msg
                
        except Exception as e:
            error_message = f"Error loading data: {str(e)}"
            self.show_message(error_message, "#FF0000")
            traceback.print_exc()
            return None, None, error_message