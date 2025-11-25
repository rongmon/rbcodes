# zgui Multi-Frame Bug Fixes - Changes Summary

## Critical Bug Fix: Frame Selection Always Used EMLINEA

### Problem
When measuring redshifts from different spectroscopic frames (SCI, SCIA, SCIB, EMLINEA, EMLINEB), all measurements were recorded under the first detected frame (EMLINEA) regardless of which frame was selected in the toolbar.

### Root Cause
`linelist_selection.py` line 307 was hardcoded to use:
```python
frame_source = self.frame_sources[0]  # Always first frame!
```

### Solution
Pass canvas reference to LineListWidget so it can read the currently selected frame directly:

---

## Files Modified

### 1. **linelist_selection.py**

#### Change 1a: Constructor (line 44)
```python
# BEFORE:
def __init__(self):

# AFTER:
def __init__(self, canvas=None):
    self.canvas = canvas  # Reference to MplCanvas
```

#### Change 1b: Frame Source Determination (lines 307-311)
```python
# BEFORE:
frame_source = self.frame_sources[0] if self.frame_sources else 'DEFAULT'

# AFTER:
if self.canvas and hasattr(self.canvas, 'current_frame'):
    frame_source = self.canvas.current_frame
else:
    frame_source = self.frame_sources[0] if self.frame_sources else 'DEFAULT'
```

#### Change 1c: Linelist 'NONE' Handling (lines 226-233)
```python
# BEFORE:
if s in 'NONE':
    self.send_linelist.emit(s)  # Sending string!

# AFTER:
if s in 'NONE':
    empty_df = pd.DataFrame(columns=['wave', 'name'])
    self.linelist = empty_df
    self.send_linelist.emit(empty_df)  # Sending DataFrame
```

#### Change 1d: Z Error Handling (lines 244-258)
Added extraction of frame_source from z_est signal (though not used with current fix, provides future flexibility):
```python
def _on_estZ_changed(self, newz):
    # newz now contains: [z_value, z_error, frame_source]
    # But we use canvas.current_frame instead
```

---

### 2. **main.py**

#### Change 2: Pass Canvas Reference (line 65)
```python
# BEFORE:
widget_z = LineListWidget()

# AFTER:
widget_z = LineListWidget(canvas=self.sc)
```

**Why**: Canvas (`self.sc`) is created before LineListWidget on line 57, so the reference is available.

---

### 3. **spec_plot.py**

#### Change 3a: Linelist Type Check (lines 1083-1086)
```python
# BEFORE:
if self.gauss_num == 0:
    self.guess_ion = GuessTransition(self.linelist, ...)  # Could crash if list

# AFTER:
# Check if linelist is available and not empty
if isinstance(self.linelist, list) or (hasattr(self.linelist, 'empty') and self.linelist.empty):
    self.send_message.emit('Error: No linelist loaded. Please select a linelist first.')
    return
```

#### Change 3b: Frame Combobox Signal Blocking (lines 561-584)
```python
# BEFORE:
self.frame_combobox.clear()
self.frame_combobox.addItem('SCI')
# ... (spurious signals during population) ...
self.frame_combobox.setCurrentIndex(...)

# AFTER:
self.frame_combobox.blockSignals(True)  # Prevent signals during population
self.frame_combobox.clear()
self.frame_combobox.addItem('SCI')
# ...
self.frame_combobox.blockSignals(False)
current_frame = self.frame_combobox.currentText()
if current_frame and current_frame != 'NONE':
    self.send_frame.emit(current_frame)  # Explicit signal after population
```

---

## Summary of Changes by Issue

### Issue 1: Frame Always EMLINEA
**Files**: linelist_selection.py, main.py
**Severity**: CRITICAL
**Fix**: Read canvas.current_frame directly instead of using frame_sources[0]

### Issue 2: Linelist TypeError
**Files**: linelist_selection.py, spec_plot.py
**Severity**: HIGH
**Fix**: Send empty DataFrame instead of string 'NONE'

### Issue 3: Frame Combobox Spurious Signals
**Files**: menu_toolbars.py
**Severity**: MEDIUM
**Fix**: Block signals during population and emit single explicit signal

---

## Testing

### Test Scripts Created
1. `test_frame_fix.py` - Tests frame combobox blocking and linelist handling
2. `test_frame_source_fix.py` - Tests canvas reference reading

### All Tests Pass ✓

---

## Lines Changed

```
linelist_selection.py: 44, 226-233, 244-258, 307-311  (4 changes)
main.py: 65                                             (1 change)
spec_plot.py: 1083-1086                                 (1 change)
menu_toolbars.py: 561-584                               (1 change, in blockSignals block)
```

**Total**: ~40 lines modified/added across 4 files

---

## Backward Compatibility

✓ All changes are backward compatible
✓ Canvas parameter is optional (defaults to None)
✓ Fallback mechanisms in place
✓ Existing CSV files load correctly
✓ No data format changes

---

## Documentation

Created 3 comprehensive guides:
1. **FRAME_SELECTION_FIX.md** - Technical explanation of this fix
2. **MULTI_FRAME_GUIDE.md** - Complete user guide and architecture
3. **CHANGES_SUMMARY.md** - This file

---

## Verification Steps

After applying these changes:

1. Load a FITS file with multiple frames (e.g., SCI, SCIA, EMLINEA)
2. Select "SCI" from frame dropdown
3. Measure a redshift
4. Click "Add to Table"
5. **Expected**: Row shows `z_SCI` column, primary_frame = 'SCI'
6. Select "EMLINEA" from frame dropdown
7. Measure same object again
8. Click "Add to Table"
9. **Expected**: Same row now has `z_EMLINEA` column, primary_frame updated
10. Right-click on row → "Set Primary Frame" → can choose which measurement to display

If all steps work, the bug is fixed!

---

## Error Messages That Might Appear (And are OK)

### "Error: No linelist loaded. Please select a linelist first."
- User right-clicked to select line before choosing a linelist
- This is expected and helpful - tells user what to do

### Frame shows 'DEFAULT'
- Canvas reference was None when button clicked
- System fell back to first detected frame
- Check that canvas is passed to LineListWidget constructor

---

## Future Improvements

If desired in future:
1. Could pass frame_source through z_est signal for additional validation
2. Could add logging to track frame selection events
3. Could add dropdown to manually select frame in LineListWidget
4. Could add keyboard shortcuts for frame switching

But current solution is minimal and robust!
