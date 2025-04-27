# Atomic Line Information Retrieval Tool

[Back to Main Page](../main_readme.md)


The `rb_setline` module provides functionality for retrieving atomic line information from various spectral line lists. It allows users to find spectral line data based on wavelength matching (exact or closest) or by species name.

## Main Features

- **Multiple Matching Methods**: Find lines by exact wavelength, closest wavelength, or species name
- **Comprehensive Line Lists**: Access to multiple spectral line catalogs for different astrophysical contexts
- **Efficient Caching**: Line lists are cached to improve performance with repeated queries
- **Robust Error Handling**: Comprehensive validation and error reporting
- **Detailed Logging**: Optional logging for troubleshooting and monitoring

## Dependencies

- NumPy
- Astropy (for ASCII file handling)
- Python standard libraries (logging, os, pathlib)

## Available Line Lists

The module supports numerous line lists for different astronomical applications:

- `atom`: Full atomic line list with oscillator strengths and damping constants
- `LLS`: Lines commonly observed in Lyman Limit Systems
- `DLA`: Lines commonly observed in Damped Lyman-alpha systems
- `LBG`: Lines commonly observed in Lyman Break Galaxies
- `Gal`: Galaxy vacuum wavelength lines
- `Eiger_Strong`: Strong lines in the Eiger catalog
- `Gal_Em`: Galaxy emission lines
- `Gal_Abs`: Galaxy absorption lines
- `Gal_long`: Extended catalog of galaxy emission and absorption lines
- `AGN`: Active Galactic Nuclei lines
- `HI_recomb`: Hydrogen recombination lines
- `HI_recomb_light`: Subset of hydrogen recombination lines

## Usage Examples

### Find the Closest Line to a Given Wavelength

```python
from rbcodes.IGM.rb_setline import rb_setline

# Find the closest line to 2796.3 Å
result = rb_setline(2796.3, 'closest')

print(f"Found: {result['name']} at {result['wave']} Å")
print(f"Oscillator strength: {result['fval']}")
```

### Match a Line by Exact Wavelength

```python
# Find an exact match for Lyman-alpha at 1215.67 Å
result = rb_setline(1215.67, 'Exact')

if len(result['wave']) > 0:
    print(f"Found exact match: {result['name']} at {result['wave']} Å")
    print(f"Oscillator strength: {result['fval']}")
    if 'gamma' in result:
        print(f"Damping constant: {result['gamma']}")
else:
    print("No exact match found")
```

### Match a Line by Species Name

```python
# Find all lines for a specific species by name
result = rb_setline(0, 'Name', target_name='HI 1215')

if len(result['wave']) > 0:
    print(f"Found {len(result['wave'])} matches for {result['name'][0]}")
    for i in range(len(result['wave'])):
        print(f"  {result['wave'][i]} Å, f = {result['fval'][i]}")
else:
    print("No matches found")
```

### Using a Different Line List

```python
# Find a line in the LBG line list
result = rb_setline(1550.0, 'closest', linelist='LBG')

print(f"Found in LBG catalog: {result['name']} at {result['wave']} Å")
```

## Return Value Structure

The function returns a dictionary with the following keys:

- `wave`: Rest frame wavelength(s) of the matched line(s) in Å
- `fval`: Oscillator strength value(s)
- `name`: Species name(s)
- `gamma`: Radiation damping constant (only available with the 'atom' line list)

For 'closest' method, results are single values. For 'Exact' and 'Name' methods, results may contain multiple values as NumPy arrays.

## API Reference

```python
rb_setline(lambda_rest, method, linelist='atom', target_name=None)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| lambda_rest | float | Rest frame wavelength in Å to match (ignored when method='Name') |
| method | str | Matching method: 'closest', 'Exact', or 'Name' |
| linelist | str | Line list to use (default: 'atom') |
| target_name | str | Required when method='Name'. Species name to match. |

## Best Practices

1. **Line Selection**:
   - Use 'closest' method for exploratory analysis to find nearby lines
   - Use 'Exact' method when you know the precise wavelength (within 0.001 Å)
   - Use 'Name' method when targeting specific transitions

2. **Working with Different Line Lists**:
   - Use 'atom' for general atomic transitions with accurate oscillator strengths
   - Use specialized lists ('LLS', 'DLA', etc.) for specific astrophysical contexts
   - Consider using `HI_recomb` for studies of Hydrogen recombination

3. **Error Handling**:
   - Always check if the result contains valid data by examining the length of the 'wave' array
   - When using 'closest' method, verify the wavelength difference is acceptable
   - When using 'Name' method, ensure consistent naming format (e.g., 'HI 1215')

4. **Performance Considerations**:
   - Line lists are automatically cached to improve performance
   - For batch processing, make repeated calls to `rb_setline()` rather than reloading line lists

## Example Workflow: Finding Multiple Lines

```python
from rbcodes.IGM.rb_setline import rb_setline
import numpy as np

# List of species we're interested in
species_list = [
    'HI 1215',    # Lyman-alpha
    'CIV 1548',   # Carbon IV doublet (first line)
    'CIV 1550',   # Carbon IV doublet (second line)
    'MgII 2796',  # Magnesium II doublet (first line)
    'MgII 2803'   # Magnesium II doublet (second line)
]

# Create a dictionary to store the results
line_info = {}

# Loop through each species and get the line information
for species in species_list:
    result = rb_setline(0, 'Name', target_name=species)
    
    if len(result['wave']) > 0:
        line_info[species] = {
            'wavelength': result['wave'][0],
            'f_value': result['fval'][0]
        }
        if 'gamma' in result:
            line_info[species]['gamma'] = result['gamma'][0]
    else:
        print(f"Warning: Could not find line data for {species}")

# Print the results
for species, info in line_info.items():
    print(f"{species}: λ = {info['wavelength']:.2f} Å, f = {info['f_value']:.4f}")
```
