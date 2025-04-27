# Spectral Rebinning Tool

The `rb_specbin` function is a specialized tool for rebinning one-dimensional spectral data to a new integer pixel scale. It handles flux, error, and wavelength arrays, providing a clean and efficient implementation for spectroscopic data processing.

## Main Features

- **Integer Pixel Rebinning**: Combines adjacent pixels into larger bins for improved signal-to-noise
- **Error Propagation**: Correctly propagates uncertainties during rebinning
- **Wavelength Handling**: Optional rebinning of associated wavelength arrays
- **Robust Input Validation**: Comprehensive error checking for input data
- **Flexible Return Format**: Dictionary-based output for easy access to rebinned data

## Dependencies

- NumPy
- Python standard math library

## Usage Examples

### Basic Rebinning of Flux Data

```python
import numpy as np
from rbcodes.IGM.rb_specbin import rb_specbin

# Create example data
flux = np.random.normal(1.0, 0.1, 1000)  # Simulated flux with noise

# Rebin by a factor of 5
result = rb_specbin(flux, 5)

# Access the rebinned flux
rebinned_flux = result['flux']

print(f"Original length: {len(flux)}")
print(f"Rebinned length: {len(rebinned_flux)}")
```

### Rebinning with Wavelength and Error Arrays

```python
import numpy as np
from rbcodes.IGM.rb_specbin import rb_specbin

# Create example data
wave = np.linspace(4000, 5000, 1000)  # Wavelength in Angstroms
flux = np.random.normal(1.0, 0.1, 1000)  # Flux with noise
error = np.ones_like(flux) * 0.1  # Constant error

# Rebin by a factor of 4
result = rb_specbin(flux, 4, var=error**2, wave=wave)

# Access the rebinned data
rebinned_wave = result['wave']
rebinned_flux = result['flux']
rebinned_error = result['error']

# Plot the original and rebinned spectra
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.errorbar(wave, flux, yerr=error, alpha=0.3, label='Original')
plt.errorbar(rebinned_wave, rebinned_flux, yerr=rebinned_error, 
             fmt='o', label='Rebinned')
plt.legend()
plt.xlabel('Wavelength (Å)')
plt.ylabel('Flux')
plt.title(f'Spectral Rebinning (nbin={4})')
plt.show()
```

### Handling Non-Divisible Array Lengths

```python
# When the array length is not divisible by nbin, the function
# automatically handles the remainder properly

# Create example data with a non-divisible length
flux = np.random.normal(1.0, 0.1, 1030)  # Length not divisible by 7

# Rebin by a factor of 7
result = rb_specbin(flux, 7)

# The result handles the remainder appropriately
print(f"Original length: {len(flux)}")
print(f"Rebinned length: {len(result['flux'])}")
print(f"Expected length: {len(flux) // 7 + 1}")
```

## Visual Example

Rebinning can significantly improve signal-to-noise ratio at the cost of spectral resolution. The diagram below illustrates this trade-off:

![Spectral Rebinning Example](/images/spectral_rebinning_example.png)

## Function Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| flux | array-like | Input flux array to be rebinned |
| nbin | int | Number of pixels to combine in each bin |
| var | array-like, optional | Input variance array (error²) |
| wave | array-like, optional | Input wavelength array |

## Return Value Structure

The function returns a dictionary with the following keys:

- `flux`: Rebinned flux array
- `error`: Rebinned error array (only if `var` was provided)
- `wave`: Rebinned wavelength array (only if `wave` was provided)

## Best Practices

1. **Choosing the Bin Size**:
   - Balance between noise reduction and preserving spectral features
   - For emission/absorption line studies, ensure bin size doesn't wash out narrow features
   - For continuum studies, larger bins can improve SNR significantly

2. **Error Handling**:
   - Provide variance (error²) rather than error directly
   - The function returns standard error (not variance) in the output dictionary

3. **Wavelength Considerations**:
   - Rebinned wavelength values are calculated as the mean of the original wavelengths
   - For non-linear wavelength scales, consider rebinning in log(wavelength) space

4. **Performance Tips**:
   - Pre-filter bad pixels before rebinning to avoid propagating poor data
   - For very large spectra, consider rebinning in sections to improve memory efficiency

5. **Validation**:
   - Always verify that spectral features of interest are preserved after rebinning
   - Calculate pre- and post-rebinning SNR to confirm improvement

## Common Use Cases

1. **Improving Signal-to-Noise Ratio**:
   - Useful for faint objects or low-exposure spectra
   - Helps reveal continuum shape in noisy data

2. **Matching Resolution**:
   - Rebin high-resolution spectra to match lower-resolution data
   - Useful for combining observations from different instruments

3. **Computational Efficiency**:
   - Reduce data volume for faster processing
   - Particularly helpful for Monte Carlo simulations with many spectral iterations

4. **Visualization Enhancement**:
   - Smoother appearance for presentation graphics
   - Clearer identification of broad spectral features

## Comparison with Other Methods

Unlike interpolation or convolution-based methods, `rb_specbin` performs simple arithmetic averaging of adjacent pixels, which:

1. Preserves the total flux
2. Correctly propagates noise statistics
3. Does not introduce spurious features
4. Maintains photon-counting statistics

This makes it particularly suitable for astronomical spectroscopy where flux conservation and proper error propagation are essential.

![Comparison of Rebinning Methods](/images/rebinning_methods_comparison.png)