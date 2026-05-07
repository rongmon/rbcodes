# rb_multispec Review TODO

Branch: `multispec-review`

Work through these one by one. Each item needs: (1) verify the issue is real, (2) fix it, (3) test.

---

## Confirmed Bugs

- [x] **BUG-1** `handle_convert_clicked` ignores return value of `convert_json_to_text()`
  - Location: `multispec.py` line 2257
  - Issue: `convert_json_to_text()` returns a 4-tuple `(line_list_path, absorbers_path, metadata_path, error)` but the result is assigned to `success` and treated as a bool — always truthy, so "Conversion Complete" shows even on failure.
  - Same pattern on line 2295 for `convert_text_to_json()`.
  - Fix: unpack the tuple and check `error is None`.

- [x] **BUG-2** Double render in `handle_redshift_submission` — redundant `plot_redshift_lines()` call
  - Location: `multispec.py` line 1242
  - Issue: `set_redshift_data()` already calls `plot_redshift_lines()` internally (line 574). Calling it again immediately draws twice.
  - Fix: remove the second call at line 1242.

- [x] **BUG-3** "None" line list submit is blocked by validation in `RedshiftInputWidget`
  - Location: `RedshiftInputWidget.py`
  - Issue: Submitting with linelist="None" shows a warning instead of clearing lines. The canvas handles "None" correctly if it receives the signal, but the widget prevents the signal from being emitted.
  - Fix: Allow "None" to pass through — the canvas `set_redshift_data()` already handles the clear case.

---

## Code Cleanliness

- [x] **CLEAN-1** Remove ~30 debug `print()` statements from `vStack.py`
  - Location: `vStack.py` — top-level import prints, `__init__`, key handler, canvas swap
  - Issue: Very noisy for end users; every key press, page turn, and layout attempt prints to console.
  - Fix: Remove or replace with `logging.debug()`.

- [x] **CLEAN-2** Remove duplicate `import matplotlib.pyplot as plt`
  - Location: `multispec.py` lines 7 and 24
  - Fix: Remove the second import.

- [x] **CLEAN-3** Remove ~300 lines of dead commented-out code
  - Location: `multispec.py` lines 1934–2231 (old `display_line_list()` implementation)
  - Fix: Delete the commented block entirely.

---

## Potential Logic Issues (needs verification)

- [ ] **LOGIC-1** `reconcile_linelists()` chaining bug in velocity clustering
  - Location: `utils.py`
  - Issue: When grouping lines by velocity threshold, each new line is compared only to the *last* entry in the current cluster, not to the cluster anchor/centroid. Lines that drift cumulatively beyond the threshold can still be merged.
  - Needs: write a test case with synthetic data to confirm or rule out.

- [~] **LOGIC-2** `IOManager` singleton state leaks between `MainWindow` instances — **DEFERRED**
  - Only affects multi-window `from_data()` usage; symptom is messages routing to wrong window.
  - `last_directories` sharing across windows is actually desirable UX.
  - Not worth fixing unless multi-window becomes a first-class use case.

- [x] **LOGIC-3** Fragile `vStack` canvas replacement — walks widget tree by attribute name
  - Location: `vStack.py` `__init__`
  - Issue: Finds `right_layout`/`main_splitter` by attribute inspection. Falls back silently to a separate `QDialog` if not found. Restoration path may have `AttributeError` risk.
  - Needs: trigger the fallback path and verify behavior.

---

## UX / Visualization Improvements (discuss before implementing)

- [ ] **UX-1** Y-axis autoscale to *visible* x-range only
  - Each panel rescales y to whatever flux is in the current zoom window (like IRAF `y` key).
  - Currently: global fixed y-range per spectrum regardless of zoom.

- [ ] **UX-2** Status bar showing cursor position in spectral coordinates
  - Show `λ_obs`, implied `z` (given current linelist), and `Δv` from nearest identified line.
  - Updates on `motion_notify_event`.

- [ ] **UX-3** Line label toggle (`L` key)
  - Hide/show all text labels without removing vertical lines.
  - Useful when many absorbers overlap.

- [ ] **UX-4** Right-click context menu on line labels
  - Options: "Remove this line", "Copy wavelength to clipboard", "Set as primary redshift".

- [ ] **UX-5** Pan by fixed velocity step (`[` / `]` keys)
  - Step by half a window width along wavelength axis (like xdisplay in IRAF).

- [ ] **UX-6** Drag-and-drop FITS file loading
  - Accept drag-and-drop onto the canvas or main window.
  - `setAcceptDrops(True)` + `dropEvent()` on `MainWindow`.

- [ ] **UX-7** Per-spectrum filename label in each subplot
  - Small text label (top-right corner) in each panel showing the source filename.

- [ ] **UX-8** `?` key help overlay or improved help dialog
  - A popup listing all keyboard shortcuts. Currently the help text goes to the message box only.

- [ ] **UX-9** Auto-restore last session on startup
  - Offer to reload the last saved JSON session (absorbers, line list, redshift).

---

## Notes

- Start with BUG-1, BUG-2, BUG-3 and CLEAN-1 through CLEAN-3 — these are low-risk and high-confidence.
- LOGIC items need test cases written first before fixing.
- UX items should be discussed and prioritized before any implementation.
- All changes go on the `multispec-review` branch.
