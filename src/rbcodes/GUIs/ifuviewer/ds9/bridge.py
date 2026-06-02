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
        except ImportError as _e:
            return False, f'no_pyds9: {_e}'
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
        Fetch region text from ds9.

        Returns a combined text: extraction regions in fk5 sky coordinates
        (for accurate WCS projection) plus annotation shapes (compass, vector,
        ruler, projection, line, segment) in image coordinates.

        Annotation shapes like compass and vector have no WCS representation in
        ds9 — they exist only in pixel space and are silently dropped when ds9
        saves in wcs/fk5 format.  We therefore do a separate image-coordinate
        save and append only the annotation lines to the sky-coord text.
        """
        if not self.available:
            return ''

        sel = 'selected ' if selected_only else ''

        # ----------------------------------------------------------------
        # Step 1: get sky-coordinate text (extraction regions)
        # ----------------------------------------------------------------
        sky_text = ''

        # Approach 1a: XPA get with explicit fk5 sky coords
        for cmd in (
            f'regions {sel}wcs fk5 degrees',
            f'regions {sel}fk5',
            f'regions {sel}wcs',
        ):
            try:
                text = self._ds9.get(cmd)
                if text and text.strip():
                    print(f"[DS9Bridge] sky regions via '{cmd}' ({len(text)} chars)")
                    sky_text = text
                    break
            except Exception:
                pass

        # Approach 1b: force wcs/fk5 then save to file
        if not sky_text:
            tmp = tempfile.mktemp(suffix='.reg')
            try:
                for set_cmd in ('regions system wcs',
                                'regions sky fk5',
                                'regions skyformat degrees'):
                    try:
                        self._ds9.set(set_cmd)
                    except Exception:
                        pass
                save_cmd = (f'regions save {tmp} selected' if selected_only
                            else f'regions save {tmp}')
                self._ds9.set(save_cmd)
                if os.path.exists(tmp):
                    with open(tmp, 'r') as fh:
                        sky_text = fh.read()
                    print(f"[DS9Bridge] sky regions via file save "
                          f"({len(sky_text)} chars)")
            except Exception as exc:
                print(f"[DS9Bridge] sky regions save failed: {exc}")
            finally:
                try:
                    os.unlink(tmp)
                except OSError:
                    pass

        # ----------------------------------------------------------------
        # Step 2: get image-coordinate text to capture annotation shapes
        # (compass, vector, ruler, projection, line, segment).
        # These shapes have no WCS equivalent and are dropped from fk5 saves.
        # ----------------------------------------------------------------
        _ANNOTATION_SHAPES = frozenset({
            'compass', 'vector', 'ruler', 'projection',
            'line', 'segment', 'text',
        })
        ann_text = ''
        tmp2 = tempfile.mktemp(suffix='.reg')
        try:
            for set_cmd in ('regions system image',):
                try:
                    self._ds9.set(set_cmd)
                except Exception:
                    pass
            save_cmd2 = (f'regions save {tmp2} selected' if selected_only
                         else f'regions save {tmp2}')
            self._ds9.set(save_cmd2)
            if os.path.exists(tmp2):
                import re as _re
                with open(tmp2, 'r') as fh:
                    raw = fh.read()
                # Keep only lines that are annotation shapes
                ann_lines = ['image']   # coordinate system header
                for line in raw.splitlines():
                    stripped = line.strip()
                    if not stripped or stripped.startswith('#'):
                        continue
                    m = _re.match(r'(\w+)\s*\(', stripped, _re.IGNORECASE)
                    if m and m.group(1).lower() in _ANNOTATION_SHAPES:
                        ann_lines.append(line)
                if len(ann_lines) > 1:
                    ann_text = '\n'.join(ann_lines)
                    print(f"[DS9Bridge] annotation shapes (image coords): "
                          f"{len(ann_lines) - 1} shape(s)")
        except Exception as exc:
            print(f"[DS9Bridge] image regions save failed: {exc}")
        finally:
            try:
                os.unlink(tmp2)
            except OSError:
                pass

        # Restore to WCS so ds9 display is not disturbed
        try:
            self._ds9.set('regions system wcs')
        except Exception:
            pass

        combined = sky_text
        if ann_text:
            combined = (combined.rstrip() + '\n' + ann_text + '\n')
        return combined

    def push_regions(self, reg_text):
        """
        Load region text directly into the current ds9 frame.

        Uses the ``regions`` XPA command which accepts inline text, so no
        temp file is needed.  Existing ds9 regions are preserved — new
        ones are added on top.

        Parameters
        ----------
        reg_text : str  — ds9 .reg format text (fk5 or image coords)

        Returns
        -------
        bool : True on success, False on failure.
        """
        if not self.available:
            return False
        try:
            self._ds9.set('regions', reg_text)
            return True
        except Exception:
            pass
        # Fallback: write temp file and load from disk
        tmp = tempfile.mktemp(suffix='.reg')
        try:
            with open(tmp, 'w') as fh:
                fh.write(reg_text)
            self._ds9.set(f'regions load {tmp}')
            return True
        except Exception as exc:
            print(f"[DS9Bridge] push_regions failed: {exc}")
            return False
        finally:
            try:
                os.unlink(tmp)
            except OSError:
                pass
