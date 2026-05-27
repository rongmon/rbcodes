"""
DS9Bridge — optional wrapper around pyds9.

All methods fail gracefully: if pyds9 is not installed or no ds9 is running,
every call is a no-op (or returns a safe empty value).  The GUI should check
``bridge.available`` before enabling ds9-specific controls.
"""
import os
import tempfile


class DS9Bridge:
    """Thin wrapper around pyds9.DS9().  Constructed lazily — connect() first."""

    def __init__(self):
        self._ds9 = None

    # ------------------------------------------------------------------
    # Connection
    # ------------------------------------------------------------------

    def connect(self):
        """
        Attempt to connect to a running ds9 instance.

        Returns
        -------
        (ok: bool, reason: str)
            ok=True on success.
            ok=False, reason='no_pyds9'  — pyds9 not installed.
            ok=False, reason='no_ds9'    — installed but no ds9 running.
        """
        try:
            import pyds9
        except ImportError:
            return False, 'no_pyds9'
        try:
            self._ds9 = pyds9.DS9()
            return True, ''
        except Exception:
            return False, 'no_ds9'

    def disconnect(self):
        self._ds9 = None

    @property
    def available(self):
        """True when a live ds9 connection exists."""
        return self._ds9 is not None

    # ------------------------------------------------------------------
    # Image sending
    # ------------------------------------------------------------------

    def send_image(self, image2d, header, frame=1):
        """
        Send a 2-D numpy array to ds9 frame *frame*.

        Writes to a temporary FITS file, sends it to ds9, then deletes the file.
        No-op if not connected.
        """
        if not self.available:
            return
        from astropy.io import fits

        hdu = fits.PrimaryHDU(image2d.astype('float32'), header=header)
        tmp = tempfile.mktemp(suffix='.fits')
        try:
            hdu.writeto(tmp, overwrite=True)
            self._ds9.set(f'frame {frame}')
            self._ds9.set(f'fits {tmp}')
        except Exception:
            pass
        finally:
            try:
                os.unlink(tmp)
            except OSError:
                pass

    # ------------------------------------------------------------------
    # WCS locking
    # ------------------------------------------------------------------

    def match_wcs(self):
        """Lock all ds9 frames to the same WCS.  No-op if not connected."""
        if not self.available:
            return
        try:
            self._ds9.set('frame match wcs')
            self._ds9.set('frame lock wcs')
        except Exception:
            pass

    # ------------------------------------------------------------------
    # Region import
    # ------------------------------------------------------------------

    def get_regions(self, selected_only=False):
        """
        Fetch region text from ds9 in fk5 sky coordinates.

        Always requests fk5 so that regions drawn on any ds9 image (HST,
        IFU, ground-based) can be correctly projected onto the IFU cube WCS.

        Returns empty string on failure.
        """
        if not self.available:
            return ''

        # Approach 1: XPA get with explicit fk5 sky coords
        sel = 'selected ' if selected_only else ''
        for cmd in (
            f'regions {sel}wcs fk5 degrees',
            f'regions {sel}fk5',
            f'regions {sel}wcs',
        ):
            try:
                text = self._ds9.get(cmd)
                if text and text.strip():
                    print(f"[DS9Bridge] regions via '{cmd}' ({len(text)} chars)")
                    return text
            except Exception:
                pass

        # Approach 2: set coordinate system to wcs/fk5, then save to temp file
        tmp = tempfile.mktemp(suffix='.reg')
        try:
            # Ask ds9 to output in fk5 sky coords before saving
            for set_cmd in ('regions system wcs',
                            'regions sky fk5',
                            'regions skyformat degrees'):
                try:
                    self._ds9.set(set_cmd)
                except Exception:
                    pass

            save_cmd = f'regions save {tmp} selected' if selected_only \
                       else f'regions save {tmp}'
            self._ds9.set(save_cmd)
            if os.path.exists(tmp):
                with open(tmp, 'r') as fh:
                    text = fh.read()
                print(f"[DS9Bridge] regions via file save ({len(text)} chars)")
                return text
        except Exception as exc:
            print(f"[DS9Bridge] regions save failed: {exc}")
        finally:
            try:
                os.unlink(tmp)
            except OSError:
                pass

        return ''
