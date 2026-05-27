"""
rb_ifuview — command-line entry point for the IFU Viewer GUI.

Usage
-----
  rb_ifuview                              # launch empty
  rb_ifuview cube.fits                    # load one file
  rb_ifuview cube1.fits cube2.fits img.fits   # load multiple files
  rb_ifuview --help                       # show this help
"""
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        prog='rb_ifuview',
        description=(
            'IFU Viewer — interactive QFitsView-style viewer for FITS '
            'datacubes and 2-D images.\n\n'
            'Supports KCWI, MUSE, and any generic 3-axis FITS cube.\n'
            'Multiple files can be loaded at once; switch between them '
            'in the sidebar.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  rb_ifuview                              # launch empty\n'
            '  rb_ifuview cube.fits                    # load one cube\n'
            '  rb_ifuview cube1.fits cube2.fits        # load two cubes\n'
            '  rb_ifuview cube.fits hst_image.fits     # cube + 2-D image\n'
        ),
    )

    parser.add_argument(
        'files',
        nargs='*',
        metavar='FILE',
        help='FITS datacube(s) or 2-D image(s) to load on startup.',
    )

    parser.add_argument(
        '--install',
        action='store_true',
        help='Print instructions for installing optional dependencies (pyds9 etc.) and exit.',
    )

    # --- future flags go here ---
    # parser.add_argument('--variance', metavar='VAR_FITS',
    #     help='Variance cube to assign to the first loaded cube.')
    # parser.add_argument('--instrument', metavar='NAME',
    #     choices=['kcwi', 'muse', 'auto'], default='auto',
    #     help='Force instrument loader (default: auto-detect from header).')

    args = parser.parse_args()

    if args.install:
        _print_install_instructions()
        return

    from rbcodes.GUIs.ifuviewer.main import launch_viewer
    launch_viewer(files=args.files if args.files else None)


def _print_install_instructions():
    print("""
rb_ifuview — optional dependency installation
=============================================

INCLUDED (installed with rbcodes by default)
  regions >= 0.5      Full ds9 region shape support (circle, box, ellipse,
                      polygon, annulus) with correct WCS handling.
                      pip install regions

OPTIONAL (install separately as needed)
  pyds9               Live ds9 bridge — send images to ds9, import regions,
                      cross-instrument region import (e.g. HST → IFU).
                      Requires SAOImageDS9 and XPA to be installed first.

  Install ds9 (macOS via Homebrew):
    brew install --cask saoimage-ds9
    brew install xpa

  Then install pyds9:
    pip install pyds9

  Verify the full chain:
    python -c "import pyds9; d = pyds9.DS9(); print(d.get('version'))"
  (ds9 must be running before running the above)

All optional features degrade gracefully — rb_ifuview works fully without
pyds9.  You will see a message box if you try to connect and it is not
available.
""")


if __name__ == '__main__':
    main()
