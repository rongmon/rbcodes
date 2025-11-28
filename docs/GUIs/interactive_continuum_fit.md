# Interactive Continuum Fitting

This tool provides an interactive interface for fitting continuum to spectral data, essential for absorption line analysis. You can choose between polynomial fitting or spline interpolation methods.

## Overview

Continuum fitting is a critical step in spectral analysis, enabling the detection and measurement of absorption or emission features. This tool provides two approaches:

- **Polynomial Fitting**: Fits a Legendre polynomial to the unmasked portions of the spectrum
- **Spline Fitting**: Creates a cubic spline that passes through user-defined anchor points

Both methods allow for masking of spectral features that should be excluded from the continuum fit, such as absorption or emission lines.

## Launching the Tool

```python
# Via Python
from rbcodes.GUIs.interactive_continuum_fit import launch_interactive_continuum_fit
launch_interactive_continuum_fit()
```

Or use it from within `launch_specgui` (rb_spec pipeline).

## Interface Layout

The interface consists of:

- **Main Plot Area**: Shows the original spectrum (top) and normalized spectrum after fitting (bottom)
- **Toolbar**: Contains standard matplotlib navigation tools plus specialized buttons
- **Sidebar**: Contains fitting controls and options organized in tabs

> **Tip**: The plot area uses standard matplotlib interactions: pan, zoom, etc. Make sure no toolbar buttons are active when adding masks or spline points.

## Basic Workflow

1. Select the fitting method: **Polynomial** or **Spline**
2. For polynomial fitting:
   - Add masks to exclude absorption or emission features
   - Set polynomial options (order, optimization)
   - Click "Fit Continuum" button
3. For spline fitting:
   - Add spline anchor points by left-clicking (for median points) or using 'b' key (for exact points)
   - Set spline options (degree, smoothing)
   - Click "Fit Continuum" button
4. Review the normalized spectrum in the bottom panel
5. Make adjustments as needed and refit
6. When satisfied, click "Accept & Return" to use the continuum fit

## Fitting Methods

### Polynomial Fitting

Uses Legendre polynomials to fit a smooth continuum to unmasked portions of the spectrum. This method is best for spectra with smooth continua and well-defined absorption or emission features.

**Options:**
- **Polynomial Order**: Higher orders can fit more complex shapes but may overfit
- **Use Weights**: When checked, uses error spectrum for weighted fitting
- **Auto Optimize**: Uses Bayesian Information Criterion to select the optimal order

**Workflow:**
1. Mask regions containing absorption or emission features
2. Set polynomial order or enable auto optimization
3. Click "Fit Continuum"

> **Tip**: For complex spectra, it's often better to use Auto Optimize with a reasonable range of orders (e.g., 1-6). This balances fitting quality with avoiding overfitting.

### Spline Fitting

Uses cubic spline interpolation between user-selected anchor points. This method gives you more direct control and works well for irregular continua or when polynomial fitting struggles.

**Options:**
- **Spline Degree**: Order of polynomial pieces (3 = cubic is standard)
- **Smoothing**: Higher values create smoother curves (0 = exact interpolation)
- **Median Window**: Size of window for calculating median flux when adding points
- **Window Type**: Units for median window (pixels or velocity/wavelength units)

**Adding Points:**
- **Left-click**: Adds a point at the median flux within the specified window around the click
- **Key 'b'**: Adds a point at the exact cursor position
- **Right-click**: Removes the closest spline point

> **Warning**: Spline fitting requires at least 3 anchor points. Place points on parts of the spectrum that represent the true continuum level, avoiding absorption or emission features.

## Working with Masks

Masks define regions to exclude from the continuum fit, typically absorption or emission features.

### Adding and Removing Masks

- **Left-click twice**: First click defines start of mask, second click defines end
- **Right-click twice**: Defines a region where masks will be removed (any mask overlapping with this region will be removed)
- **Auto-Mask button**: Automatically detects potential features to mask based on statistical deviations
- **Manual Entry button**: Opens a dialog to precisely specify mask boundaries
- **Reset Masks button**: Clears all masks
- **Undo Last button**: Removes the most recently added mask

> **Note**: Overlapping masks are automatically merged into a single continuous mask region.

### Auto-Masking Parameters (Advanced tab)

- **Auto-Mask Sigma**: Threshold for detecting features (in standard deviations)
- **Min Mask Width**: Minimum width for a detected feature to be masked

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `f` | Fit continuum |
| `r` | Reset masks |
| `R` (capital) | Reset everything (masks, spline points, and fit) |
| `z` | Undo last mask |
| `a` | Accept results |
| `escape` | Cancel/close without saving |
| `b` | Add spline point at exact cursor position (spline mode only) |
| `c` | Clear all spline points (spline mode only) |
| `p` | Switch to polynomial fitting mode |
| `s` | Switch to spline fitting mode |
| `m` | Open manual mask entry dialog |
| `+` / `=` | Zoom in |
| `-` / `_` | Zoom out |
| `0` | Reset zoom |
| `h` | Show help dialog |

## Tips and Best Practices

- **Fitting complex spectra**: Try both methods and compare results
- **For polynomial fitting**: Start with low orders (2-3) and increase only if needed
- **For spline fitting**: Place anchor points at regular intervals in regions representing the true continuum
- **When masking**: Be generous with mask boundaries to fully exclude features
- **Auto-optimizing**: Use a reasonable range (e.g., 1-6) to avoid overfitting
- **Check normalization**: The bottom panel shows the normalized spectrum - the continuum should follow 1.0 in regions without features
- **Use keyboard shortcuts**: They speed up the workflow significantly

## Common Issues to Avoid

- Using polynomial orders that are too high, causing overfitting
- Placing spline points on absorption or emission features instead of the true continuum
- Not masking the full extent of spectral features
- Not including enough spline points to adequately define the continuum shape

## Important

> The quality of the continuum fit directly affects the accuracy of equivalent width measurements and column density calculations. Take time to ensure the fit is appropriate for your science goals.
