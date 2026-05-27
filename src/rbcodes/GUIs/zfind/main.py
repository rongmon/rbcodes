"""
main.py — CLI entry point for rb_zfind.

Usage
-----
rb_zfind                          # empty dialog, load spectrum from within
rb_zfind spectrum.fits            # load file and open dialog
rb_zfind spectrum.fits -l ISM     # with default linelist
rb_zfind spectrum.fits -m absorption  # absorption mode
"""

import sys
import argparse


def launch_zfind(spec_or_file=None, linelist=None, mode='emission'):
    """
    Launch the ZFindDialog programmatically.

    Parameters
    ----------
    spec_or_file : str, rb_spectrum, tuple, or None
        - str       → loaded via rb_spectrum.from_file()
        - rb_spectrum → used directly
        - tuple     → passed to rb_spectrum.from_tuple()
        - None      → empty dialog (user loads from within)
    linelist : str or None
        Default linelist name to pre-select (e.g. 'ISM')
    mode : 'emission' or 'absorption'
    """
    from PyQt5.QtWidgets import QApplication
    from rbcodes.GUIs.zfind.io import to_rb_spectrum
    from rbcodes.GUIs.zfind.dialog import ZFindDialog

    app = QApplication.instance() or QApplication(sys.argv[:1])

    spec = None
    if spec_or_file is not None:
        spec = to_rb_spectrum(spec_or_file)

    dialog = ZFindDialog(spec=spec, mode=mode, default_linelist=linelist)
    dialog.exec_()


def main():
    """Entry point for rb_zfind console script."""
    parser = argparse.ArgumentParser(
        description='rb_zfind: Semi-automated redshift and absorber finder',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
rb_zfind                           # open empty dialog
rb_zfind spectrum.fits             # load and open
rb_zfind spectrum.fits -l ISM      # with default linelist
rb_zfind spectrum.fits -m absorption   # absorption mode
"""
    )
    parser.add_argument('fitsfile', nargs='?', default=None,
                        help='FITS file to analyze (optional)')
    parser.add_argument('-l', '--linelist', default=None,
                        help='Default linelist to pre-select (e.g. ISM)')
    parser.add_argument('-m', '--mode', default='emission',
                        choices=['emission', 'absorption'],
                        help='Search mode (default: emission)')
    args = parser.parse_args()

    launch_zfind(args.fitsfile, linelist=args.linelist, mode=args.mode)


if __name__ == '__main__':
    main()
