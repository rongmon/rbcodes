"""
Region parsing, rasterization, and overlay drawing for the IFU Viewer.

Primary path: uses the ``regions`` package (pip install regions) which supports
all ds9 format shapes with correct WCS handling and matplotlib integration.

Fallback path (no ``regions`` installed): hand-rolled numpy rasterizer that
produces correct masks but draws only a cross marker on the image.
"""
import re
import numpy as np
from matplotlib.path import Path

# ---------------------------------------------------------------------------
# Optional dependency
# ---------------------------------------------------------------------------

try:
    from regions import Regions as _Regions, PixelRegion as _PixelRegion, \
        SkyRegion as _SkyRegion
    _HAS_REGIONS = True
except ImportError:
    _HAS_REGIONS = False

# Shapes that produce extraction masks
EXTRACTION_SHAPES = {'circle', 'box', 'ellipse', 'polygon', 'annulus'}
# Shapes that are annotation only (no extraction mask)
ANNOTATION_SHAPES = {'text', 'compass', 'ruler', 'projection', 'vector',
                     'point', 'line', 'segment'}
# Annotation shapes that the ``regions`` package does NOT parse — must always
# be picked up by the regex fallback even when the regions package is available.
_UNHANDLED_BY_REGIONS = frozenset({'compass', 'ruler', 'projection', 'vector',
                                    'line', 'segment'})


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_ds9_regions(reg_text):
    """
    Parse ds9 region format text → list of region dicts.

    Each dict:
    {
      'shape'    : str              e.g. 'circle'
      'args'     : list[str]        raw args (used by fallback only)
      'system'   : str              'fk5', 'icrs', 'image', …
      'name'     : str              from ``# text={…}`` or ''
      '_reg_obj' : Region or None   native regions object (preferred path)
    }

    Excluded regions (leading ``-``) are silently dropped.
    Annotation shapes (compass, vector, ruler, etc.) are included with
    extract=False so they can be displayed as overlays.
    """
    if _HAS_REGIONS:
        result = _parse_with_regions(reg_text)
        if result is None:
            # Parse hard-failed → full fallback
            return _parse_fallback(reg_text)
        # The regions package silently drops compass, vector, ruler, projection,
        # line, segment.  Always supplement with the regex fallback for those.
        extra = [r for r in _parse_fallback(reg_text)
                 if r.get('shape', '') in _UNHANDLED_BY_REGIONS]
        if extra:
            print(f"[spatial_mask] supplemented with {len(extra)} "
                  f"annotation shape(s) not handled by regions package")
        return result + extra
    return _parse_fallback(reg_text)


def region_to_mask(region, wcs, ny, nx):
    """
    Convert one region dict → boolean mask (ny, nx).

    Uses the ``regions`` package when available (accurate for all shapes).
    Falls back to hand-rolled numpy for the common shapes.
    Clips silently to image bounds; returns all-False if region is outside FOV.
    """
    reg_obj = region.get('_reg_obj')
    if _HAS_REGIONS and reg_obj is not None:
        try:
            return _mask_from_reg_obj(reg_obj, wcs, ny, nx)
        except Exception as exc:
            print(f"[spatial_mask] _mask_from_reg_obj failed: {exc}")
    # fallback
    return _mask_fallback(region, wcs, ny, nx)


def draw_region_overlay(region, ax, wcs, color=None, draw_label=True):
    """
    Draw the region shape on a matplotlib axes.

    Uses ``regions`` → correct shape (circle, ellipse, box, polygon, annulus).
    Falls back to a cross marker when the package is not available.

    Parameters
    ----------
    region     : dict — from parse_ds9_regions()
    ax         : matplotlib Axes
    wcs        : astropy WCS (2D) or None
    color      : str or None — override color; uses file color if None
    draw_label : bool — draw text label from region.meta['text']

    Returns
    -------
    list of matplotlib artists added to *ax*
    """
    artists = []

    reg_obj = region.get('_reg_obj')
    if _HAS_REGIONS and reg_obj is not None:
        try:
            artists.extend(_draw_with_regions(reg_obj, ax, wcs, color))
        except Exception:
            pass

    if not artists:
        artists.extend(_draw_fallback_marker(region, ax, wcs, color))

    # Text label
    if draw_label:
        name = region.get('name', '')
        if name:
            try:
                artists.extend(_draw_label(region, ax, wcs, color, ny=None, nx=None))
            except Exception:
                pass

    return artists


def region_center_sky(region, wcs):
    """
    Return (ra_deg, dec_deg) of the region centre.

    Returns (None, None) if conversion fails or no WCS available.
    """
    reg_obj = region.get('_reg_obj')

    # Preferred: extract from native object
    if _HAS_REGIONS and reg_obj is not None:
        try:
            if isinstance(reg_obj, _SkyRegion):
                ctr = reg_obj.center
                return float(ctr.ra.deg), float(ctr.dec.deg)
            elif isinstance(reg_obj, _PixelRegion) and wcs is not None:
                ctr = reg_obj.center
                sky = wcs.pixel_to_world(ctr.x, ctr.y)
                return float(sky.ra.deg), float(sky.dec.deg)
        except Exception:
            pass

    # Fallback: parse args manually
    if wcs is None:
        return None, None
    try:
        args   = region.get('args', [])
        system = region.get('system', 'fk5')
        is_sky = system not in ('image', 'physical', 'linear')
        if is_sky:
            return float(args[0]), float(args[1])
        cx = float(args[0]) - 1.0
        cy = float(args[1]) - 1.0
        ra, dec = wcs.all_pix2world([[cx, cy]], 0)[0]
        return float(ra), float(dec)
    except Exception:
        return None, None


def iau_name(ra_deg, dec_deg, prefix='spec_'):
    """
    Build an IAU-format filename from sky coordinates.

    Example: iau_name(152.12, +21.56, 'spec_') → 'spec_J100828.8+213336.fits'
    Falls back to decimal coords if astropy.coordinates unavailable.
    """
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        s = c.to_string('hmsdms', sep='', precision=1).replace(' ', '')
        return f'{prefix}{s}.fits'
    except Exception:
        return f'{prefix}RA{ra_deg:.4f}_Dec{dec_deg:+.4f}.fits'


def _sep_arcsec(wcs, x0, y0, x1, y1):
    """
    Great-circle angular separation in arcsec between pixels (x0,y0) and (x1,y1).

    Calls all_pix2world with origin=0 (0-indexed).  For celestial axes astropy
    always returns degrees, so the haversine formula below is correct regardless
    of the CUNIT header values.  Returns None on failure.
    """
    try:
        ra0, dec0 = wcs.all_pix2world([[x0, y0]], 0)[0]
        ra1, dec1 = wcs.all_pix2world([[x1, y1]], 0)[0]
        d_ra  = np.deg2rad(ra1 - ra0) * np.cos(np.deg2rad((dec0 + dec1) / 2.0))
        d_dec = np.deg2rad(dec1 - dec0)
        return np.rad2deg(np.sqrt(d_ra ** 2 + d_dec ** 2)) * 3600.0
    except Exception:
        return None


def _radius_arcsec(wcs, cx, cy, r_pix):
    """
    Convert a pixel radius to arcsec by measuring 1-pixel steps along x and y
    and scaling.  Averaging x and y handles IFUs with different plate scales on
    each axis (e.g. KCWI slices vs. along-slice direction).
    """
    sx = _sep_arcsec(wcs, cx, cy, cx + r_pix, cy)
    sy = _sep_arcsec(wcs, cx, cy, cx,          cy + r_pix)
    vals = [v for v in (sx, sy) if v is not None]
    return sum(vals) / len(vals) if vals else None


def _pixel_axis_pa(wcs, cx, cy, dx=0, dy=1):
    """
    Position angle (degrees East of North, CCW, 0–360) of a pixel-space
    direction, computed via astropy SkyCoord.position_angle.

    Using SkyCoord.position_angle is essential — a manual -(ΔRA)*cos(dec)
    formula has the wrong sign when CDELT1 > 0 (RA increasing East, as some
    IFU pipelines produce), giving the mirror angle (e.g. 40.5° instead of
    the correct 319.5°).

    dx, dy define the direction in pixel space:
      dx=0, dy=1  →  pixel y-axis (height direction) — default for ds9 box
      dx=1, dy=0  →  pixel x-axis (width  direction)

    Returns 0.0 on failure.
    """
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        ra0, dec0 = wcs.all_pix2world([[cx,      cy     ]], 0)[0]
        ra1, dec1 = wcs.all_pix2world([[cx + dx, cy + dy]], 0)[0]
        sky0 = SkyCoord(ra=ra0 * u.deg, dec=dec0 * u.deg)
        sky1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg)
        return float(sky0.position_angle(sky1).deg)   # 0–360°, E of N, CCW
    except Exception:
        return 0.0


def regions_to_ds9_text(aperture_list, wcs=None):
    """
    Convert a list of aperture dicts → ds9 .reg format text.

    Each dict:
      type    : 'single' | 'circle' | 'annulus' | 'rect'
      cx, cy  : pixel centre (0-indexed)
      radius  : pixels (circle/annulus)
      bg_inner, bg_outer : pixels (annulus, optional)
      x1, y1, x2, y2 : pixels (rect, 0-indexed)
      label   : str

    Sizes are converted to arcsec via direct WCS pixel-to-sky calls so that
    the result is correct regardless of CDELT units or non-square IFU plate
    scales.
    """
    lines = [
        '# Region file format: DS9 version 4.1',
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" '
        'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1',
    ]
    lines.append('fk5' if wcs is not None else 'image')

    for ap in aperture_list:
        atype   = ap.get('type', 'single')
        label   = ap.get('label', '')
        comment = f' # text={{{label}}}' if label else ''

        # Pixel centre (0-indexed)
        cx = float(ap.get('cx', (ap.get('x1', 0) + ap.get('x2', 0)) / 2.0))
        cy = float(ap.get('cy', (ap.get('y1', 0) + ap.get('y2', 0)) / 2.0))

        ra = dec = None
        if wcs is not None:
            try:
                ra, dec = wcs.all_pix2world([[cx, cy]], 0)[0]
            except Exception:
                pass

        if atype == 'single':
            if ra is not None:
                lines.append(f'circle({ra:.6f},{dec:.6f},1.0"){comment}')
            else:
                lines.append(f'circle({cx + 1:.1f},{cy + 1:.1f},1){comment}')

        elif atype == 'circle':
            r      = float(ap.get('radius', 3))
            bg_in  = ap.get('bg_inner')
            bg_out = ap.get('bg_outer')
            if ra is not None:
                r_as = _radius_arcsec(wcs, cx, cy, r)
                if r_as is not None:
                    lines.append(
                        f'circle({ra:.6f},{dec:.6f},{r_as:.3f}"){comment}')
                    if bg_in is not None and bg_out is not None:
                        ri_as = _radius_arcsec(wcs, cx, cy, float(bg_in))
                        ro_as = _radius_arcsec(wcs, cx, cy, float(bg_out))
                        if ri_as is not None and ro_as is not None:
                            lines.append(
                                f'annulus({ra:.6f},{dec:.6f},'
                                f'{ri_as:.3f}",{ro_as:.3f}"){comment}')
            else:
                lines.append(
                    f'circle({cx + 1:.1f},{cy + 1:.1f},{r:.1f}){comment}')
                if bg_in is not None and bg_out is not None:
                    lines.append(
                        f'annulus({cx + 1:.1f},{cy + 1:.1f},'
                        f'{float(bg_in):.1f},{float(bg_out):.1f}){comment}')

        elif atype == 'annulus':
            r_in  = float(ap.get('bg_inner', ap.get('radius', 3)))
            r_out = float(ap.get('bg_outer', r_in + 3))
            if ra is not None:
                ri_as = _radius_arcsec(wcs, cx, cy, r_in)
                ro_as = _radius_arcsec(wcs, cx, cy, r_out)
                if ri_as is not None and ro_as is not None:
                    lines.append(
                        f'annulus({ra:.6f},{dec:.6f},'
                        f'{ri_as:.3f}",{ro_as:.3f}"){comment}')
            else:
                lines.append(
                    f'annulus({cx + 1:.1f},{cy + 1:.1f},'
                    f'{r_in:.1f},{r_out:.1f}){comment}')

        elif atype == 'rect':
            x1 = float(ap.get('x1', cx))
            y1 = float(ap.get('y1', cy))
            x2 = float(ap.get('x2', cx))
            y2 = float(ap.get('y2', cy))
            if ra is not None:
                # Width: sky distance along the pixel x-axis at constant y = cy
                # Height: sky distance along the pixel y-axis at constant x = cx
                w_as = _sep_arcsec(wcs, x1, cy, x2, cy)
                h_as = _sep_arcsec(wcs, cx, y1, cx, y2)
                # PA of pixel y-axis (height direction) East of North, CCW.
                # ds9 WCS box angle=0 means height→North; angle rotates the
                # whole box CCW from that default orientation, so we need the
                # angle of the *height* (pixel y) axis from North, not width.
                pa   = _pixel_axis_pa(wcs, cx, cy, dx=0, dy=1)
                if w_as is not None and h_as is not None:
                    lines.append(
                        f'box({ra:.6f},{dec:.6f},'
                        f'{w_as:.3f}",{h_as:.3f}",{pa:.4f}){comment}')
            else:
                # ds9 image box: centre (1-indexed) + full width + full height
                cxds9 = (x1 + x2) / 2.0 + 1
                cyds9 = (y1 + y2) / 2.0 + 1
                w = abs(x2 - x1)
                h = abs(y2 - y1)
                lines.append(
                    f'box({cxds9:.1f},{cyds9:.1f},{w:.1f},{h:.1f},0){comment}')

    return '\n'.join(lines) + '\n'


# ---------------------------------------------------------------------------
# regions-package implementation
# ---------------------------------------------------------------------------

def _parse_with_regions(reg_text):
    """Parse using astropy-regions; return list of dicts or None on failure."""
    import warnings
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            reg_list = _Regions.parse(reg_text, format='ds9')
    except Exception as exc:
        print(f"[spatial_mask] regions.parse failed: {exc}")
        return None

    result = []
    for reg in reg_list:
        shape = _reg_shape_name(reg)
        extract = bool(shape)  # annotation types return '' → display only
        if not shape:
            # Unknown/annotation type — include as display-only overlay
            # Use the class name as the shape label for drawing fallback
            shape = type(reg).__name__.lower().replace('region', '').strip('_') or 'unknown'
            print(f"[spatial_mask] annotation/unknown region type: {type(reg).__name__} → display only")

        name = ''
        meta = getattr(reg, 'meta', {})
        name = str(meta.get('text', '') or meta.get('label', '') or '')

        system = 'image' if isinstance(reg, _PixelRegion) else 'fk5'

        result.append({
            'shape':    shape,
            'args':     [],      # not needed when _reg_obj present
            'system':   system,
            'name':     name,
            '_reg_obj': reg,
            'extract':  extract,
        })

    n_ext = sum(1 for r in result if r.get('extract', True))
    n_ann = len(result) - n_ext
    print(f"[spatial_mask] regions package found {n_ext} extraction region(s), {n_ann} annotation(s)")
    return result


def _reg_shape_name(reg):
    """Map a regions class → shape name string, or '' if annotation/unsupported."""
    cls = type(reg).__name__.lower()
    # Longer prefixes must be checked before shorter ones (e.g. circleannulus before circle)
    _MAP = [
        ('circleannulus',    'annulus'),
        ('ellipseannulus',   'annulus'),
        ('rectangleannulus', 'annulus'),
        ('circle',           'circle'),
        ('ellipse',          'ellipse'),
        ('rectangle',        'box'),
        ('polygon',          'polygon'),
        ('point',            ''),    # annotation
        ('text',             ''),
    ]
    for prefix, shape in _MAP:
        if cls.startswith(prefix):
            return shape
    return ''


def _mask_from_reg_obj(reg_obj, wcs, ny, nx):
    """Rasterize a regions object → boolean (ny, nx) mask."""
    if isinstance(reg_obj, _SkyRegion):
        if wcs is None:
            return np.zeros((ny, nx), dtype=bool)
        pix_reg = reg_obj.to_pixel(wcs)
    else:
        pix_reg = reg_obj

    mask_obj = pix_reg.to_mask(mode='center')
    img = mask_obj.to_image((ny, nx))
    if img is None:
        ctr = getattr(pix_reg, 'center', None)
        print(f"[spatial_mask] to_image returned None — pixel center={ctr}, image size=({ny},{nx})")
        return np.zeros((ny, nx), dtype=bool)
    return np.asarray(img, dtype=bool)


def _draw_with_regions(reg_obj, ax, wcs, color):
    """Draw reg_obj on ax; return list of artists added."""
    if isinstance(reg_obj, _SkyRegion):
        if wcs is None:
            return []
        pix_reg = reg_obj.to_pixel(wcs)
    else:
        pix_reg = reg_obj

    # Get file color if caller didn't override
    vis = getattr(reg_obj, 'visual', {})
    file_color = vis.get('edgecolor') or vis.get('facecolor')
    draw_color = color or file_color or '#a6e3a1'

    artist = pix_reg.as_artist()
    artist.set_facecolor('none')
    artist.set_edgecolor(draw_color)
    artist.set_linewidth(1.4)
    artist.set_zorder(11)
    artist.set_picker(5)
    ax.add_artist(artist)
    return [artist]


def _draw_label(region, ax, wcs, color, ny, nx):
    """Draw text label at region centroid. Returns list of text artists."""
    name = region.get('name', '')
    if not name:
        return []

    reg_obj = region.get('_reg_obj')
    cx = cy = None

    if _HAS_REGIONS and reg_obj is not None:
        try:
            if isinstance(reg_obj, _SkyRegion) and wcs is not None:
                pix = reg_obj.to_pixel(wcs)
                ctr = pix.center
                cx, cy = float(ctr.x), float(ctr.y)
            elif isinstance(reg_obj, _PixelRegion):
                ctr = reg_obj.center
                cx, cy = float(ctr.x), float(ctr.y)
        except Exception:
            pass

    if cx is None:
        # Fallback: centroid from args
        args = region.get('args', [])
        try:
            cx = float(args[0]) - (0 if region.get('system') in ('fk5', 'icrs') else 1)
            cy = float(args[1]) - (0 if region.get('system') in ('fk5', 'icrs') else 1)
        except Exception:
            return []

    txt = ax.text(cx, cy + 1, name,
                  color=color or '#a6e3a1', fontsize=7,
                  va='bottom', ha='center', zorder=13,
                  bbox=dict(facecolor='#1e1e2e', alpha=0.5, pad=1, edgecolor='none'))
    return [txt]


# ---------------------------------------------------------------------------
# Fallback: hand-rolled numpy rasterizer + matplotlib patch drawing
# ---------------------------------------------------------------------------

def _parse_fallback(reg_text):
    """Parse ds9 region text without the regions package → list of dicts."""
    print("[spatial_mask] using fallback regex parser")
    regions_out = []
    coord_system = 'fk5'

    for raw in reg_text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('#'):
            # ds9 writes annotation shapes as comment lines, e.g.:
            #   # compass(ra,dec,r") compass=icrs {N} {E} 1 1
            #   # vector(x,y,len,angle) vector
            # Strip the leading '#' and keep only if it's a known annotation.
            rest = line[1:].strip()
            m_ann = re.match(r'(\w+)\s*\(', rest, re.IGNORECASE)
            if m_ann and m_ann.group(1).lower() in ANNOTATION_SHAPES:
                line = rest   # treat as a regular shape line below
            else:
                continue      # ordinary comment — discard

        comment = ''
        if '#' in line:
            line, comment = line.split('#', 1)
            line    = line.strip()
            comment = comment.strip()

        lower = line.lower()
        if lower in ('fk5', 'fk4', 'icrs', 'galactic', 'ecliptic',
                     'image', 'physical', 'wcs', 'linear'):
            coord_system = lower
            continue
        if lower.startswith('global ') or lower == 'global':
            continue

        excluded = line.startswith('-')
        if excluded:
            line = line[1:].strip()

        m = re.match(r'(\w+)\s*\(([^)]+)\)', line, re.IGNORECASE)
        if not m:
            continue

        shape = m.group(1).lower()
        args  = [a.strip() for a in m.group(2).split(',')]

        if excluded:
            continue
        if shape not in EXTRACTION_SHAPES and shape not in ANNOTATION_SHAPES:
            continue  # truly unknown shape

        name = ''
        nm   = re.search(r'text\s*=\s*\{([^}]*)\}', comment, re.IGNORECASE)
        if nm:
            name = nm.group(1).strip()

        extractable = shape in EXTRACTION_SHAPES
        regions_out.append({
            'shape':    shape,
            'args':     args,
            'system':   coord_system,
            'name':     name,
            '_reg_obj': None,
            'extract':  extractable,
        })

    n_ext = sum(1 for r in regions_out if r.get('extract', True))
    n_ann = len(regions_out) - n_ext
    print(f"[spatial_mask] fallback found {n_ext} extraction region(s), {n_ann} annotation(s)")
    return regions_out


def _mask_fallback(region, wcs, ny, nx):
    """Rasterize a region dict using numpy (no regions package required)."""
    shape  = region.get('shape', '')
    args   = region.get('args', [])
    system = region.get('system', 'fk5')
    is_sky = system not in ('image', 'physical', 'linear')

    try:
        if shape == 'circle':
            return _mask_circle(args, wcs, ny, nx, is_sky)
        if shape == 'box':
            return _mask_box(args, wcs, ny, nx, is_sky)
        if shape == 'ellipse':
            return _mask_ellipse(args, wcs, ny, nx, is_sky)
        if shape == 'polygon':
            return _mask_polygon(args, wcs, ny, nx, is_sky)
        if shape == 'annulus':
            return _mask_annulus(args, wcs, ny, nx, is_sky)
    except Exception:
        pass
    return np.zeros((ny, nx), dtype=bool)


def _draw_fallback_marker(region, ax, wcs, color):
    """
    Draw region shape using matplotlib patches directly from args.

    Handles circle, ellipse, box, polygon, annulus in both pixel (physical/image)
    and sky coordinates.  Falls back to a cross if the shape is not recognised.
    """
    from matplotlib.patches import Circle, Ellipse, Rectangle, Polygon as MplPolygon
    from matplotlib.transforms import Affine2D

    draw_color = color or '#a6e3a1'
    shape  = region.get('shape', '')
    args   = region.get('args', [])
    system = region.get('system', 'fk5')
    is_sky = system not in ('image', 'physical', 'linear')
    artists = []

    def _to_pix(ra_or_x, dec_or_y):
        """Convert to pixel coords (0-indexed)."""
        if is_sky and wcs is not None:
            xy = wcs.all_world2pix([[float(ra_or_x), float(dec_or_y)]], 0)
            return float(xy[0][0]), float(xy[0][1])
        return float(ra_or_x) - 1.0, float(dec_or_y) - 1.0

    def _sz(s):
        """Size arg → pixels."""
        if is_sky and wcs is not None:
            return _parse_sky_size(s) / _pix_scale_deg(wcs)
        return _strip_unit(s)

    patch = None
    try:
        if shape == 'circle' and len(args) >= 3:
            cx, cy = _to_pix(args[0], args[1])
            r = _sz(args[2])
            patch = Circle((cx, cy), r, fill=False,
                           edgecolor=draw_color, lw=1.4, zorder=11, picker=5)

        elif shape == 'ellipse' and len(args) >= 4:
            cx, cy = _to_pix(args[0], args[1])
            r1 = _sz(args[2])
            r2 = _sz(args[3])
            angle = float(args[4]) if len(args) > 4 else 0.0
            # Ellipse(xy, width, height) — full axes, not semi-axes
            patch = Ellipse((cx, cy), 2 * r1, 2 * r2,
                            angle=angle, fill=False,
                            edgecolor=draw_color, lw=1.4, zorder=11, picker=5)

        elif shape == 'box' and len(args) >= 4:
            cx, cy = _to_pix(args[0], args[1])
            w = _sz(args[2])
            h = _sz(args[3])
            angle = float(args[4]) if len(args) > 4 else 0.0
            # Rectangle anchor is bottom-left; rotate around centre
            rect = Rectangle((-w / 2, -h / 2), w, h, fill=False,
                             edgecolor=draw_color, lw=1.4, zorder=11, picker=5)
            t = (Affine2D().rotate_deg(angle)
                           .translate(cx, cy)
                 + ax.transData)
            rect.set_transform(t)
            ax.add_patch(rect)
            artists.append(rect)

        elif shape == 'polygon' and len(args) >= 6:
            coords = [float(a) for a in args]
            xs, ys = coords[0::2], coords[1::2]
            if is_sky and wcs is not None:
                pix = wcs.all_world2pix(list(zip(xs, ys)), 0)
                verts = [(p[0], p[1]) for p in pix]
            else:
                verts = [(x - 1, y - 1) for x, y in zip(xs, ys)]
            patch = MplPolygon(verts, closed=True, fill=False,
                               edgecolor=draw_color, lw=1.4, zorder=11, picker=5)

        elif shape == 'annulus' and len(args) >= 4:
            cx, cy = _to_pix(args[0], args[1])
            r_in  = _sz(args[2])
            r_out = _sz(args[3])
            for r in (r_in, r_out):
                ring = Circle((cx, cy), r, fill=False,
                              edgecolor=draw_color, lw=1.4, zorder=11, picker=5)
                ax.add_patch(ring)
                artists.append(ring)

        # ---- annotation shapes ----

        elif shape in ('ruler', 'line', 'segment', 'projection') and len(args) >= 4:
            # Two-endpoint shapes: (x1,y1,x2,y2[,width])
            x1p, y1p = _to_pix(args[0], args[1])
            x2p, y2p = _to_pix(args[2], args[3])
            ln, = ax.plot([x1p, x2p], [y1p, y2p], '-',
                          color=draw_color, lw=1.4, zorder=11, picker=5)
            artists.append(ln)
            # endpoints
            for px, py in ((x1p, y1p), (x2p, y2p)):
                mk, = ax.plot(px, py, 'o', color=draw_color, ms=3,
                              zorder=12, picker=5)
                artists.append(mk)

        elif shape == 'vector' and len(args) >= 3:
            # vector(x,y,length,angle) — angle CCW from +x in image coords
            from matplotlib.patches import FancyArrowPatch
            cx, cy = _to_pix(args[0], args[1])
            length = _sz(args[2])
            angle  = float(args[3]) if len(args) > 3 else 0.0
            tip_x  = cx + length * np.cos(np.deg2rad(angle))
            tip_y  = cy + length * np.sin(np.deg2rad(angle))
            arr = FancyArrowPatch((cx, cy), (tip_x, tip_y),
                                  arrowstyle='->', color=draw_color,
                                  linewidth=1.4, mutation_scale=10, zorder=11)
            ax.add_patch(arr)
            artists.append(arr)

        elif shape == 'compass' and len(args) >= 2:
            # compass(ra,dec,length) — draw North and East arms from centre.
            # ds9 saves compass as a comment line:
            #   # compass(ra,dec,r") compass=icrs {N} {E} 1 1
            # so is_sky=True when parsed from an fk5 region file.
            from matplotlib.patches import FancyArrowPatch

            cx, cy = _to_pix(args[0], args[1])

            if wcs is not None and is_sky and len(args) >= 3:
                # _parse_sky_size('20.675"') returns degrees → * 3600 = arcsec
                arm_as = _parse_sky_size(args[2]) * 3600.0
            elif len(args) >= 3:
                arm_as = _strip_unit(args[2])   # image coords: already pixels
            else:
                arm_as = 15.0

            # Compute N and E tip pixel positions via SkyCoord offsets
            tips = {}  # 'N' → (tip_x, tip_y), 'E' → (tip_x, tip_y)
            if wcs is not None and is_sky:
                try:
                    from astropy.coordinates import SkyCoord
                    import astropy.units as u
                    ra0, dec0 = wcs.all_pix2world([[cx, cy]], 0)[0]
                    cen_sky = SkyCoord(ra=ra0 * u.deg, dec=dec0 * u.deg)
                    for pa_deg, lbl in ((0.0, 'N'), (90.0, 'E')):
                        tip_sky = cen_sky.directional_offset_by(
                            pa_deg * u.deg, arm_as * u.arcsec)
                        tip_pix = wcs.all_world2pix(
                            [[tip_sky.ra.deg, tip_sky.dec.deg]], 0)[0]
                        tx, ty = float(tip_pix[0]), float(tip_pix[1])
                        if np.isfinite(tx) and np.isfinite(ty):
                            tips[lbl] = (tx, ty)
                except Exception as e:
                    print(f"[spatial_mask] compass tip calc failed: {e}")

            # Fallback for any arm that failed WCS projection
            if 'N' not in tips:
                tips['N'] = (cx, cy + arm_as)
            if 'E' not in tips:
                tips['E'] = (cx + arm_as, cy)

            mk, = ax.plot(cx, cy, 'o', color=draw_color, ms=3, zorder=12)
            artists.append(mk)
            for lbl, (tip_x, tip_y) in tips.items():
                arr = FancyArrowPatch((cx, cy), (tip_x, tip_y),
                                     arrowstyle='->', color=draw_color,
                                     linewidth=1.4, mutation_scale=10,
                                     zorder=11)
                ax.add_patch(arr)
                artists.append(arr)
                txt = ax.text(tip_x, tip_y, lbl, color=draw_color,
                              fontsize=7, ha='center', va='bottom', zorder=13)
                artists.append(txt)

    except Exception as exc:
        print(f"[spatial_mask] fallback draw failed ({shape}): {exc}")

    if patch is not None:
        ax.add_patch(patch)
        artists.append(patch)

    # If nothing was drawn, fall back to a cross at the region centre
    if not artists:
        try:
            # For regions-package objects with no args, try to get centre from _reg_obj
            if not args:
                reg_obj = region.get('_reg_obj')
                if reg_obj is not None and _HAS_REGIONS:
                    pix = (reg_obj.to_pixel(wcs)
                           if isinstance(reg_obj, _SkyRegion) and wcs else reg_obj)
                    ctr = getattr(pix, 'center', None)
                    if ctr is not None:
                        args = [ctr.x, ctr.y]
            if not args:
                return artists
            cx, cy = _to_pix(args[0], args[1])
            mk, = ax.plot(cx, cy, '+', color=draw_color, ms=10, mew=1.6,
                          zorder=12, picker=5)
            artists.append(mk)
        except Exception:
            pass

    return artists


# ---------------------------------------------------------------------------
# Numpy rasterization helpers (fallback only)
# ---------------------------------------------------------------------------

def _sky_to_pix(ra, dec, wcs):
    xy = wcs.all_world2pix([[ra, dec]], 0)
    return float(xy[0][0]), float(xy[0][1])


def _pix_scale_deg(wcs):
    from astropy.wcs.utils import proj_plane_pixel_scales
    return float(np.mean(np.abs(proj_plane_pixel_scales(wcs))))


def _parse_sky_size(s):
    s = s.strip()
    if s.endswith('"'):  return float(s[:-1]) / 3600.0
    if s.endswith("'"): return float(s[:-1]) / 60.0
    if s.lower().endswith('d'): return float(s[:-1])
    if s.lower().endswith('r'): return float(s[:-1]) * 180.0 / np.pi
    return float(s)


def _strip_unit(s):
    return float(re.sub(r'[^\d.eE+\-]', '', s.strip()))


def _mask_circle(args, wcs, ny, nx, is_sky):
    if is_sky and wcs is not None:
        cx, cy = _sky_to_pix(float(args[0]), float(args[1]), wcs)
        r_pix  = _parse_sky_size(args[2]) / _pix_scale_deg(wcs)
    else:
        cx, cy, r_pix = float(args[0])-1, float(args[1])-1, _strip_unit(args[2])
    yy, xx = np.mgrid[0:ny, 0:nx]
    return (xx-cx)**2 + (yy-cy)**2 <= r_pix**2


def _mask_box(args, wcs, ny, nx, is_sky):
    if is_sky and wcs is not None:
        cx, cy = _sky_to_pix(float(args[0]), float(args[1]), wcs)
        psc    = _pix_scale_deg(wcs)
        w_pix  = _parse_sky_size(args[2]) / psc
        h_pix  = _parse_sky_size(args[3]) / psc
    else:
        cx, cy = float(args[0])-1, float(args[1])-1
        w_pix, h_pix = _strip_unit(args[2]), _strip_unit(args[3])
    angle_rad = np.deg2rad(_strip_unit(args[4]) if len(args) > 4 else 0)
    ca, sa = np.cos(angle_rad), np.sin(angle_rad)
    yy, xx = np.mgrid[0:ny, 0:nx]
    dx, dy = xx-cx, yy-cy
    xr =  dx*ca + dy*sa
    yr = -dx*sa + dy*ca
    return (np.abs(xr) <= w_pix/2) & (np.abs(yr) <= h_pix/2)


def _mask_ellipse(args, wcs, ny, nx, is_sky):
    if is_sky and wcs is not None:
        cx, cy = _sky_to_pix(float(args[0]), float(args[1]), wcs)
        psc    = _pix_scale_deg(wcs)
        r1, r2 = _parse_sky_size(args[2])/psc, _parse_sky_size(args[3])/psc
    else:
        cx, cy = float(args[0])-1, float(args[1])-1
        r1, r2 = _strip_unit(args[2]), _strip_unit(args[3])
    angle_rad = np.deg2rad(_strip_unit(args[4]) if len(args) > 4 else 0)
    ca, sa = np.cos(angle_rad), np.sin(angle_rad)
    yy, xx = np.mgrid[0:ny, 0:nx]
    dx, dy = xx-cx, yy-cy
    xr =  dx*ca + dy*sa
    yr = -dx*sa + dy*ca
    return (xr/max(r1,1e-9))**2 + (yr/max(r2,1e-9))**2 <= 1.0


def _mask_polygon(args, wcs, ny, nx, is_sky):
    coords = [float(a) for a in args]
    xs, ys = coords[0::2], coords[1::2]
    if is_sky and wcs is not None:
        pix = wcs.all_world2pix(list(zip(xs, ys)), 0)
        pxs, pys = [p[0] for p in pix], [p[1] for p in pix]
    else:
        pxs, pys = [x-1 for x in xs], [y-1 for y in ys]
    verts = list(zip(pxs, pys)) + [(pxs[0], pys[0])]
    path  = Path(verts, closed=True)
    yy, xx = np.mgrid[0:ny, 0:nx]
    pts = np.column_stack([xx.ravel().astype(float), yy.ravel().astype(float)])
    return path.contains_points(pts).reshape(ny, nx)


def _mask_annulus(args, wcs, ny, nx, is_sky):
    if is_sky and wcs is not None:
        cx, cy = _sky_to_pix(float(args[0]), float(args[1]), wcs)
        psc    = _pix_scale_deg(wcs)
        r_in   = _parse_sky_size(args[2]) / psc
        r_out  = _parse_sky_size(args[3]) / psc
    else:
        cx, cy  = float(args[0])-1, float(args[1])-1
        r_in, r_out = _strip_unit(args[2]), _strip_unit(args[3])
    yy, xx = np.mgrid[0:ny, 0:nx]
    d2 = (xx-cx)**2 + (yy-cy)**2
    return (d2 >= r_in**2) & (d2 <= r_out**2)
