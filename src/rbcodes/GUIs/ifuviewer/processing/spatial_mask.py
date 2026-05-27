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
# Shapes that are annotation only — drop silently
ANNOTATION_SHAPES = {'text', 'compass', 'ruler', 'projection', 'vector',
                     'point', 'line', 'segment'}


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

    Excluded regions (leading ``-``) and annotation shapes are silently dropped.
    """
    if _HAS_REGIONS:
        result = _parse_with_regions(reg_text)
        # Only skip fallback if regions package found at least one region
        if result:
            return result
        # result is None (parse error) or [] (no extraction shapes matched)
        # → always fall through to fallback so we surface something useful
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


def regions_to_ds9_text(aperture_list, wcs=None):
    """
    Convert a list of aperture dicts → ds9 .reg format text.

    Each dict:
      type    : 'single' | 'circle' | 'annulus' | 'rect'
      cx, cy  : pixel centre (0-indexed)
      radius  : pixels (circle/annulus)
      bg_inner, bg_outer : pixels (annulus, optional)
      x1, y1, x2, y2 : pixels (rect)
      label   : str
    """
    lines = [
        '# Region file format: DS9 version 4.1',
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" '
        'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1',
    ]
    lines.append('fk5' if wcs is not None else 'image')

    for ap in aperture_list:
        atype  = ap.get('type', 'single')
        label  = ap.get('label', '')
        comment = f' # text={{{label}}}' if label else ''
        ra = dec = pscale = None

        if wcs is not None:
            try:
                cx = ap.get('cx', (ap.get('x1', 0) + ap.get('x2', 0)) / 2)
                cy = ap.get('cy', (ap.get('y1', 0) + ap.get('y2', 0)) / 2)
                ra, dec = wcs.all_pix2world([[cx, cy]], 0)[0]
                from astropy.wcs.utils import proj_plane_pixel_scales
                pscale = float(np.mean(np.abs(proj_plane_pixel_scales(wcs)))) * 3600
            except Exception:
                pass

        if atype == 'single':
            if ra is not None:
                lines.append(f'circle({ra:.6f},{dec:.6f},1.0"){comment}')
            else:
                cx, cy = ap['cx'] + 1, ap['cy'] + 1
                lines.append(f'circle({cx},{cy},1){comment}')

        elif atype == 'circle':
            r = ap.get('radius', 3)
            if ra is not None:
                lines.append(f'circle({ra:.6f},{dec:.6f},{r * pscale:.2f}"){comment}')
            else:
                cx, cy = ap['cx'] + 1, ap['cy'] + 1
                lines.append(f'circle({cx},{cy},{r}){comment}')

        elif atype == 'annulus':
            r_in  = ap.get('bg_inner', ap.get('radius', 3))
            r_out = ap.get('bg_outer', r_in + 3)
            if ra is not None:
                lines.append(
                    f'annulus({ra:.6f},{dec:.6f},'
                    f'{r_in * pscale:.2f}",{r_out * pscale:.2f}"){comment}')
            else:
                cx, cy = ap['cx'] + 1, ap['cy'] + 1
                lines.append(f'annulus({cx},{cy},{r_in},{r_out}){comment}')

        elif atype == 'rect':
            if ra is not None:
                try:
                    ra2, dec2 = wcs.all_pix2world([[ap['x2'], ap['y2']]], 0)[0]
                    w = abs(ra2 - ra) * 3600 * np.cos(np.deg2rad(dec))
                    h = abs(dec2 - dec) * 3600
                    lines.append(
                        f'box({ra:.6f},{dec:.6f},{w:.2f}",{h:.2f}",0){comment}')
                except Exception:
                    cx = (ap['x1'] + ap['x2']) / 2 + 1
                    cy = (ap['y1'] + ap['y2']) / 2 + 1
                    w  = ap['x2'] - ap['x1']
                    h  = ap['y2'] - ap['y1']
                    lines.append(f'box({cx},{cy},{w},{h},0){comment}')
            else:
                cx = (ap['x1'] + ap['x2']) / 2 + 1
                cy = (ap['y1'] + ap['y2']) / 2 + 1
                w  = ap['x2'] - ap['x1']
                h  = ap['y2'] - ap['y1']
                lines.append(f'box({cx},{cy},{w},{h},0){comment}')

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
        if not shape:
            # Log unknown shapes so we can extend _MAP if needed
            print(f"[spatial_mask] skipping region type: {type(reg).__name__}")
            continue

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
        })

    print(f"[spatial_mask] regions package found {len(result)} extraction region(s)")
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
        if not line or line.startswith('#'):
            continue

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

        if excluded or shape in ANNOTATION_SHAPES or shape not in EXTRACTION_SHAPES:
            continue

        name = ''
        nm   = re.search(r'text\s*=\s*\{([^}]*)\}', comment, re.IGNORECASE)
        if nm:
            name = nm.group(1).strip()

        regions_out.append({
            'shape':    shape,
            'args':     args,
            'system':   coord_system,
            'name':     name,
            '_reg_obj': None,
        })

    print(f"[spatial_mask] fallback found {len(regions_out)} extraction region(s)")
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

    except Exception as exc:
        print(f"[spatial_mask] fallback draw failed ({shape}): {exc}")

    if patch is not None:
        ax.add_patch(patch)
        artists.append(patch)

    # If nothing was drawn, fall back to a cross
    if not artists:
        try:
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
