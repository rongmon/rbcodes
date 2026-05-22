"""
test_ds9.py — Verify pyds9 + ds9 communication chain.

Run this WHILE ds9 is already open.

Usage:
    python test_ds9.py

Requirements for success:
    - SAOImageDS9 installed and on PATH (or XPA_METHOD=local set)
    - pyds9 installed  (pip install pyds9)
    - ds9 window open before running this script
"""

import sys


def check(label, fn):
    """Run fn(), print PASS/FAIL with label. Return result or None on failure."""
    try:
        result = fn()
        print(f"  [PASS] {label}")
        return result
    except Exception as exc:
        print(f"  [FAIL] {label}")
        print(f"         {type(exc).__name__}: {exc}")
        return None


def main():
    print("=" * 60)
    print("ds9 / pyds9 verification script")
    print("=" * 60)

    # ------------------------------------------------------------------ #
    # 1. Import pyds9
    # ------------------------------------------------------------------ #
    print("\n[1] Import pyds9")
    pyds9 = check("import pyds9", lambda: __import__("pyds9"))
    if pyds9 is None:
        print("\n  HINT: pip install pyds9")
        print("        Also make sure XPA is installed (brew install xpa on macOS).")
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # 2. Connect to running ds9
    # ------------------------------------------------------------------ #
    print("\n[2] Connect to running ds9")

    def _connect():
        d = pyds9.DS9()
        return d

    d = check("pyds9.DS9() — connect to running ds9", _connect)
    if d is None:
        print("\n  HINTS:")
        print("    - Make sure ds9 is open BEFORE running this script.")
        print("    - On macOS, ds9 may not be on PATH. Try:")
        print("        export PATH=$PATH:/Applications/SAOImageDS9.app/Contents/MacOS")
        print("    - XPA issues? Try: export XPA_METHOD=local")
        print("      then restart both ds9 and this script.")
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # 3. Get ds9 version
    # ------------------------------------------------------------------ #
    print("\n[3] Query ds9 version")
    version = check("d.get('version')", lambda: d.get("version").strip())
    if version:
        print(f"         ds9 version: {version}")

    # ------------------------------------------------------------------ #
    # 4. Send a test image
    # ------------------------------------------------------------------ #
    print("\n[4] Send a synthetic 100×100 noise image to ds9 frame 1")

    def _send_image():
        import numpy as np

        rng = np.random.default_rng(42)
        data = rng.standard_normal((100, 100)).astype("float32")

        d.set("frame 1")
        # set_np2arr sends the array directly over XPA — no temp file needed
        d.set_np2arr(data)
        d.set("zoom to fit")
        d.set("scale zscale")

    check("d.set_np2arr() — send array directly to ds9", _send_image)

    # ------------------------------------------------------------------ #
    # 5. Round-trip: read back the image size from ds9
    # ------------------------------------------------------------------ #
    print("\n[5] Round-trip — read image dimensions back from ds9")

    def _roundtrip():
        width  = d.get("fits width").strip()
        height = d.get("fits height").strip()
        assert width == "100" and height == "100", (
            f"Expected 100×100, got {width}×{height}"
        )
        return width, height

    dims = check("d.get('fits width/height') == 100×100", _roundtrip)
    if dims:
        print(f"         Reported size: {dims[0]} × {dims[1]}")

    # ------------------------------------------------------------------ #
    # 6. Verify XPA can list ds9 processes
    # ------------------------------------------------------------------ #
    print("\n[6] List XPA access points (xpaaccess)")
    import subprocess

    def _xpa_list():
        result = subprocess.run(
            ["xpaaccess", "-n", "ds9"],
            capture_output=True, text=True, timeout=5
        )
        count = result.stdout.strip()
        assert count != "0", "xpaaccess reports 0 ds9 instances"
        return count

    count = check("xpaaccess -n ds9", _xpa_list)
    if count:
        print(f"         ds9 XPA access points found: {count}")

    # ------------------------------------------------------------------ #
    # Done
    # ------------------------------------------------------------------ #
    print("\n" + "=" * 60)
    print("All checks complete.  If all steps show [PASS], ds9 + pyds9 are")
    print("working and Phase 10 of the IFU viewer is safe to implement.")
    print("=" * 60)


if __name__ == "__main__":
    main()
