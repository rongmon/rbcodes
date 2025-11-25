# Frame Selection Bug Fix - Summary

## Problem Statement

When measuring redshifts from different spectroscopic frames (SCI, SCIA, SCIB, EMLINEA, EMLINEB, etc.), all measurements were being recorded under a single frame column (always EMLINEA as primary frame), regardless of which frame was actually selected in the toolbar.

**User Report**: "It loads different frames, I can toggle between frames and do different measurements. It always shows EMLINEA as primary no matter what frame is initiated."

## Root Cause

The bug was in **linelist_selection.py** at line 307:

```python
# OLD CODE (BROKEN)
frame_source = self.frame_sources[0] if self.frame_sources else 'DEFAULT'
```

This line **always used the FIRST frame** from the detected list (which happened to be EMLINEA), regardless of which frame was currently selected in the toolbar.

### Why This Happened

The signal flow was:
1. User selects frame in toolbar (e.g., 'SCI')
2. Canvas updates: `self.current_frame = 'SCI'`
3. User fits Gaussian or enters z value
4. Canvas emits: `send_z_est.emit([z, z_err, 'SCI'])`
5. LineListWidget receives but **doesn't use the frame from the signal**
6. When "Add to Table" button clicked, LineListWidget has **forgotten which frame was active**
7. Falls back to: `self.frame_sources[0]` → Always EMLINEA

## Solution

Instead of trying to pass the frame through the signal chain, **LineListWidget now reads the current frame directly from the canvas** at the moment the button is clicked:

### Changes Made

#### 1. **linelist_selection.py:44** - Updated `__init__` signature
```python
# NEW CODE
def __init__(self, canvas=None):
    super().__init__()
    self.canvas = canvas  # Reference to MplCanvas to read current_frame
    # ... rest of initialization
```

#### 2. **linelist_selection.py:307-311** - Updated frame source determination
```python
# NEW CODE
if self.canvas and hasattr(self.canvas, 'current_frame'):
    frame_source = self.canvas.current_frame
else:
    frame_source = self.frame_sources[0] if self.frame_sources else 'DEFAULT'
```

#### 3. **main.py:65** - Pass canvas reference when creating widget
```python
# OLD CODE
widget_z = LineListWidget()

# NEW CODE
widget_z = LineListWidget(canvas=self.sc)
```

## How It Works Now

1. User selects frame in toolbar (e.g., 'SCI')
2. Canvas updates: `self.current_frame = 'SCI'`
3. User measures redshift (manual entry, Gaussian fit, etc.)
4. User clicks "Add to Table" button
5. **LineListWidget reads directly**: `self.canvas.current_frame` → 'SCI'
6. Data sent to table with: `frame_source: 'SCI'`
7. Table creates column: `z_SCI` and `z_err_SCI`
8. Primary frame set to: 'SCI'

## Benefits of This Approach

✓ **Simplicity**: Reads value directly from source of truth (canvas) at moment of use
✓ **No Signal Chain Issues**: No need to pass frame through multiple signals
✓ **Thread Safe**: Gets current value at the right time, not a cached value
✓ **Handles Frame Changes**: If user changes frame before clicking button, the NEW frame is used
✓ **Fallback Support**: If canvas reference missing, falls back to first available frame

## Testing

Test script: `test_frame_source_fix.py`

Validates:
- ✓ Widget reads 'SCI' when canvas.current_frame = 'SCI'
- ✓ Widget reads 'EMLINEA' when canvas.current_frame = 'EMLINEA'
- ✓ Widget reads 'SCIA' when canvas.current_frame = 'SCIA'
- ✓ Fallback works when canvas is None
- ✓ Fallback works when canvas lacks current_frame attribute

## Expected User Experience Now

1. Load FITS file with frames: SCI, SCIA, SCIB, EMLINEA, EMLINEB, EMLINE, CONT
2. Select frame "SCI" from toolbar dropdown
3. Measure redshift using any method
4. Click "Add to Table"
5. **Result**: Entry created with `z_SCI` column, primary_frame = 'SCI'
6. Select frame "EMLINEA" from toolbar dropdown
7. Measure redshift for same object
8. Click "Add to Table"
9. **Result**: Same row updated with `z_EMLINEA` column, primary_frame updated to 'EMLINEA'
10. Right-click on row in table → "Set Primary Frame" → can select which measurement to show as primary

## Files Modified

1. **linelist_selection.py** - Added canvas parameter and use it to read current_frame
2. **main.py** - Pass canvas reference when creating LineListWidget
3. **spec_plot.py** - Already had proper frame selection logic and signal emission (no changes needed)
4. **tableview_pandas.py** - Already properly handles frame_source in data dict (no changes needed)

## Files NOT Modified (But Contain Related Logic)

- **menu_toolbars.py** - Frame combobox signal blocking (fixed in earlier commit)
- **gui_io.py** - Frame detection logic (working correctly)
- **spec_plot.py** - Canvas frame tracking (already correct)

## Edge Cases Handled

- No canvas reference: Falls back to `frame_sources[0]`
- Canvas exists but no `current_frame` attribute: Falls back to `frame_sources[0]`
- Empty frame_sources list: Uses 'DEFAULT'
- Frame changed between z measurement and button click: Uses the current selection

## Backwards Compatibility

✓ All existing code continues to work
✓ CSV files with existing data load correctly
✓ No changes to data format or signal signatures
✓ Canvas parameter is optional (defaults to None)
