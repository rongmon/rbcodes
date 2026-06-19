"""
main.py — CLI entry point for rb_zfind.

Usage
-----
rb_zfind                          # empty dialog, load spectrum from within
rb_zfind spectrum.fits            # load file and open dialog
rb_zfind spectrum.fits -l zfind_galaxy  # with default linelist
"""

import sys
import os
import argparse
os.environ['QT_AUTO_SCREEN_SCALE_FACTOR'] = '0'
os.environ['QT_ENABLE_HIGHDPI_SCALING'] = '0'


def launch_zfind(spec_or_file=None, linelist=None, **kwargs):
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
        Default linelist to pre-select (e.g. 'zfind_em', 'zfind_galaxy')
    """
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtCore import Qt
    from rbcodes.GUIs.zfind.io import to_rb_spectrum
    from rbcodes.GUIs.zfind.dialog import ZFindDialog

    QApplication.setAttribute(Qt.AA_DisableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication.instance() or QApplication(sys.argv[:1])

    spec = None
    if spec_or_file is not None:
        spec = to_rb_spectrum(spec_or_file)

    dialog = ZFindDialog(spec=spec, default_linelist=linelist)
    dialog.exec_()


def main():
    """Entry point for rb_zfind console script."""
    parser = argparse.ArgumentParser(
        description='rb_zfind: Semi-automated redshift finder',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
rb_zfind                                  # open empty dialog
rb_zfind spectrum.fits                    # load and open
rb_zfind spectrum.fits -l zfind_galaxy    # with default linelist
"""
    )
    parser.add_argument('fitsfile', nargs='?', default=None,
                        help='FITS file to analyze (optional)')
    parser.add_argument('-l', '--linelist', default=None,
                        help='Default linelist (e.g. zfind_em, zfind_galaxy)')
    args = parser.parse_args()

    launch_zfind(args.fitsfile, linelist=args.linelist)


if __name__ == '__main__':
    main()
