# Multi-Frame Redshift Measurement System - Complete Guide

## Overview

The zgui now supports measuring redshifts from multiple spectroscopic frames within the same FITS file and tracking which frame's measurement is "primary" in the results table.

### What is a "Frame"?

In spectroscopy, a FITS file can contain multiple exposures or measurements:
- **SCI**: Main science spectrum
- **SCIA**, **SCIB**: Science variants or different exposures
- **EMLINEA**, **EMLINEB**: Emission line centered exposures
- **EMLINE**: Generic emission line exposure
- **CONT**: Continuum spectrum
- And others...

## User Workflow

### Step 1: Load a FITS File with Multiple Frames

```bash
rb_zgui -f myfile.fits
```

The system automatically detects all available frames in the FITS file.

### Step 2: Select a Frame

Use the **frame dropdown** in the toolbar (top right area):
```
┌─────────────┐
│ SCI ▼       │  ← Dropdown shows: SCI, SCIA, SCIB, EMLINEA, EMLINEB, etc.
└─────────────┘
```

**Default**: SCI (if available)

### Step 3: Measure Redshift in Selected Frame

Choose measurement method:

#### Method A: Manual Entry
1. Enter estimated redshift in "Estimated z" field
2. Click "Add to Table"
3. **Result**: Measurement recorded in `z_SCI` (or selected frame)

#### Method B: Gaussian Fit
1. Adjust Gaussian parameters
2. Right-click on spectrum to select emission line
3. Redshift auto-calculated
4. Click "Add to Table"
5. **Result**: Measurement recorded in `z_SCI` with source = 'Gaussian'

#### Method C: Multi-Gaussian Fit
1. Set `#Gauss` to 2 or more
2. Fit multiple Gaussians
3. Click "Add to Table"
4. **Result**: Measurement recorded in `z_SCI` with source = 'Multi-Gaussian'

### Step 4: Switch Frames and Measure Again

1. Select different frame from dropdown (e.g., "EMLINEA")
2. Measure redshift using any method
3. Click "Add to Table"
4. **Result**: Same object row gets `z_EMLINEA` column, same row updated

### Step 5: Set Primary Frame

Right-click on row in table:
```
┌────────────────────────────┐
│ Set Primary Frame          │
├────────────────────────────┤
│ ✓ SCI        (current)     │
│   SCIA                     │
│   EMLINEA                  │
└────────────────────────────┘
```

**What this does**:
- The selected frame's z/z_err values become the main "z" and "z_err" in the table
- A ★ indicator shows which frame is primary
- CSV export includes all frame measurements

### Step 6: Export Results

Button: **Save**

CSV file contains:
```
Name,z,z_err,z_source,primary_frame,z_SCI,z_err_SCI,z_SCIA,z_err_SCIA,z_EMLINEA,z_err_EMLINEA,...
obj1,0.501,0.005,Gaussian,EMLINEA,0.500,0.010,0.502,0.008,0.501,0.005
obj2,0.801,0.010,Manual,SCI,0.801,0.010,,,0.805,0.015
```

## Architecture

### Key Components

#### 1. Frame Detection (gui_io.py)

```python
def _detect_frame_sources(filepath):
    """Detect all available frames in FITS file"""
    # Looks at both HDU names and header keywords
    # Returns: ['SCI', 'SCIA', 'SCIB', 'EMLINEA', 'EMLINEB', 'EMLINE', 'CONT']
```

**Priorities**:
1. HDU names (most reliable): FITS HDU names exactly match frame names
2. Header keywords: EMLINEA_*, EMLINEB_* patterns
3. Ignores generic names: SCI, FLUX, ERROR (already handled)

#### 2. Frame Combobox (menu_toolbars.py:561-584)

```python
def _select_file(self, i):
    # ... file loading ...

    if self.mW.toggle_frames:
        # Populate frame dropdown
        self.frame_combobox.blockSignals(True)  # Prevent spurious signals
        self.frame_combobox.clear()

        # Add SCI first (main frame)
        if 'SCI' in self.frames:
            self.frame_combobox.addItem('SCI')

        # Add all others alphabetically
        for frame_name in sorted(self.frames.keys()):
            if frame_name != 'SCI':
                self.frame_combobox.addItem(frame_name)

        self.frame_combobox.blockSignals(False)
        # Emit signal with selected frame
        self.send_frame.emit(self.frame_combobox.currentText())
```

**Why blockSignals?** During population, each `addItem()` triggers `currentTextChanged`. We block these intermediate signals and emit ONE explicit signal with the final frame.

#### 3. Canvas Frame Tracking (spec_plot.py:58, 1187-1192)

```python
# Initialization
self.current_frame = 'SCI'  # Track active frame

# Signal receiver
def _on_frame_selected(self, frame_name):
    """Called when user selects frame from toolbar"""
    self.current_frame = frame_name
    print(f"Canvas now displaying: {frame_name}")
```

#### 4. Frame Source Determination (linelist_selection.py:307-311)

**THE KEY FIX**: When "Add to Table" button clicked, read the current frame directly from canvas:

```python
def _on_button_clicked(self, sfilename):
    # ... prepare data ...

    # Read current frame directly from canvas (source of truth)
    if self.canvas and hasattr(self.canvas, 'current_frame'):
        frame_source = self.canvas.current_frame  # e.g., 'SCI'
    else:
        frame_source = self.frame_sources[0] if self.frame_sources else 'DEFAULT'

    data = {
        'Name': filename,
        'z': measured_z,
        'z_err': measured_err,
        'frame_source': frame_source,  # This is the key!
        'z_source': 'Manual'/'Gaussian'/'Multi-Gaussian',
        ...
    }
    self.send_data.emit(data)
```

**Why this works**:
- LineListWidget doesn't need to remember the frame
- Just reads what the canvas CURRENTLY shows
- Gets the most up-to-date value at the moment of button click
- If user changed frames before clicking, new frame is used

#### 5. Table Data Storage (tableview_pandas.py:72-117)

```python
def _on_sent_data(self, sent_data):
    """Receive measurement data and store in table"""

    # Extract frame name from data
    frame_source = sent_data.pop('frame_source', 'DEFAULT')

    # Create frame-specific columns if needed
    self._ensure_frame_columns(frame_source)

    # Update or create row
    if row_exists:
        # Update existing row
        s['z'] = sent_data['z']
        s['z_err'] = sent_data['z_err']
        s['z_SCI'] = sent_data['z']  # Store in frame column too
        s['z_err_SCI'] = sent_data['z_err']
        s['primary_frame'] = frame_source  # Most recent is primary
    else:
        # Create new row
        sent_data['z_SCI'] = sent_data['z']
        sent_data['z_err_SCI'] = sent_data['z_err']
        sent_data['primary_frame'] = frame_source
        self.estZ = self.estZ.append(sent_data)
```

#### 6. Primary Frame Selection (tableview_pandas.py:165-242)

```python
def _on_context_menu(self, position):
    """Right-click context menu to change primary frame"""

    # Find available frames for this row (non-NaN z values)
    available_frames = []
    for frame_name in self.frame_columns.keys():
        z_col = f'z_{frame_name}'
        if pd.notna(row_data[z_col]):
            available_frames.append(frame_name)

    # Create menu with checkmark for current primary
    for frame_name in available_frames:
        action = menu.addAction(frame_name)
        if frame_name == row_data['primary_frame']:
            action.setText(f"✓ {frame_name}")  # Mark current
        action.triggered.connect(lambda f=frame_name: self._set_primary_frame(row, f))

def _set_primary_frame(self, row, frame_name):
    """Update z/z_err columns with values from selected frame"""

    # Copy frame-specific values to main z/z_err columns
    z_col = f'z_{frame_name}'
    z_err_col = f'z_err_{frame_name}'
    self.estZ.iloc[row]['z'] = self.estZ.iloc[row][z_col]
    self.estZ.iloc[row]['z_err'] = self.estZ.iloc[row][z_err_col]
    self.estZ.iloc[row]['primary_frame'] = frame_name
```

## Signal Flow Diagram

```
┌──────────────────────────────────────────────────────────────┐
│                    MEASUREMENT WORKFLOW                      │
└──────────────────────────────────────────────────────────────┘

   User selects frame "SCI"
            ↓
      Toolbar.frame_combobox
            ↓
      send_frame.emit('SCI')
            ↓
      Canvas._on_frame_selected('SCI')
            ↓
      self.current_frame = 'SCI'
            ↓
      User measures z (Gaussian, Manual, etc.)
            ↓
      Canvas.send_z_est.emit([z_value, z_error])
            ↓
      LineListWidget._on_estZ_changed()
            ↓
      User clicks "Add to Table"
            ↓
      LineListWidget._on_button_clicked()
            ↓
      frame_source = self.canvas.current_frame  ← READS DIRECTLY FROM CANVAS
            ↓
      data['frame_source'] = 'SCI'
            ↓
      send_data.emit(data)
            ↓
      Table._on_sent_data(data)
            ↓
      Creates z_SCI column
      Sets primary_frame = 'SCI'
            ↓
      Table displays with ★ indicator
```

## Column Structure in Results Table

### Standard Columns (always present)
- `Name` - Object identifier
- `z` - Redshift (value from primary frame)
- `z_err` - Error (from primary frame)
- `z_source` - How z was determined: 'Manual'/'Gaussian'/'Multi-Gaussian'
- `Confidence` - User confidence rating
- `Linelist` - Which line list was used
- `Flag` - User comments
- `primary_frame` - Which frame is primary (shown with ★)

### Dynamic Frame Columns (created per frame)
- `z_SCI`, `z_err_SCI` - Measurements from SCI frame
- `z_SCIA`, `z_err_SCIA` - Measurements from SCIA frame
- `z_EMLINEA`, `z_err_EMLINEA` - Measurements from EMLINEA frame
- etc. for each detected frame

### Example Table

```
Name       z     z_err  z_source     primary_frame  z_SCI  z_err_SCI  z_SCIA  z_err_SCIA  z_EMLINEA  z_err_EMLINEA
obj1      0.501  0.005  Gaussian    ★ EMLINEA      0.500  0.010      0.502   0.008       0.501      0.005
obj2      0.801  0.010  Manual      ★ SCI          0.801  0.010      NaN     NaN         0.805      0.015
obj3      0.302  0.008  Multi-G     ★ SCIA         0.300  0.012      0.302   0.008       NaN        NaN
```

## Import/Export (CSV)

### Loading CSV with Frame Measurements

CSV file:
```csv
Name,z,z_err,z_source,primary_frame,z_SCI,z_err_SCI,z_SCIA,z_err_SCIA
obj1,0.501,0.005,Gaussian,EMLINEA,0.500,0.010,0.502,0.008
```

When loaded:
1. Table reads all columns
2. Frame columns detected automatically (z_*, z_err_*)
3. `primary_frame` value used to determine displayed z/z_err
4. All frame measurements preserved for context menu

### Saving CSV with Frame Measurements

All columns exported including:
- Main z/z_err (from primary frame)
- Frame-specific z/z_err columns
- primary_frame tracking

## Troubleshooting

### Issue: All measurements going to EMLINEA

**Cause**: Canvas not properly updated with selected frame
**Check**:
1. Does frame dropdown show different frames? (If not, frame detection failed)
2. Does dropdown respond to clicks? (If not, signal blocking issue)
3. Is frame name displayed in toolbar after selection?

**Fix**: Ensure `_on_frame_selected()` is being called in Canvas

### Issue: Can't change primary frame with context menu

**Cause**: Insufficient frame measurements (only one frame measured)
**Check**: Right-click only shows menu if 2+ frames have measurements
**Fix**: Measure in multiple frames before trying to set primary

### Issue: Frame columns not created in table

**Cause**: frame_source not in data dict when "Add to Table" clicked
**Check**: Print statement should show frame name
**Fix**: Verify canvas reference passed correctly to LineListWidget

## Key Concepts

### Immutability of Measurements
- Once a measurement is made, it's stored in frame-specific column
- Primary frame can change, but underlying measurement persists
- CSV export captures all measurements

### Primary Frame
- Most recent measurement becomes primary automatically
- User can override with right-click context menu
- Affects what is shown in main z/z_err columns
- Exported to CSV for reproducibility

### Signal Flow Simplification
- Frame information doesn't pass through signal chain
- LineListWidget reads current value from Canvas when needed
- Eliminates timing issues and signal ordering problems
- More robust to frame changes

## Testing Checklist

- [ ] Load FITS with multiple frames
- [ ] Frame dropdown populated correctly
- [ ] Can toggle between frames
- [ ] Each frame click updates display
- [ ] Measure in SCI, save to table
- [ ] Measure in EMLINEA, save to same row
- [ ] Table shows z_SCI and z_EMLINEA columns
- [ ] Primary frame indicator appears
- [ ] Right-click shows context menu
- [ ] Can change primary frame
- [ ] CSV export includes all frames
- [ ] CSV re-import loads all data correctly
- [ ] Display updates when primary frame changed

## Related Issues Fixed

### Issue 1: Frame Combobox Signal Blocking (menu_toolbars.py)
**Problem**: Multiple signals during combobox population caused spurious frame changes
**Solution**: Block signals during population, emit single explicit signal

### Issue 2: Linelist Empty DataFrame (linelist_selection.py)
**Problem**: Sending string 'NONE' instead of DataFrame caused TypeError
**Solution**: Create empty DataFrame when no linelist selected

### Issue 3: Frame Source Determination (linelist_selection.py)
**Problem**: Using `frame_sources[0]` always selected first frame
**Solution**: Read `canvas.current_frame` directly at moment of use
