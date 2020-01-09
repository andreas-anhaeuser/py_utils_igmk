"""Helpers to coordinate reference systems (CRS)."""

# standard modules
from pyproj import Proj

# local modules
from utils_igmk import units as unit_utils

def lonlat_to_utm(lon, lat, crs=None, units='m'):
    """Return x, y."""
    if crs is None:
        crs = get_utm_proj4string(lon, lat, units)

    # create projection
    if isinstance(crs, Proj):
        proj = crs
    else:
        proj = Proj(crs)

    # apply projection
    x, y = proj(lon, lat)

    # scale
    units = retrieve_units_from_crs(crs)
    scale_factor = unit_utils.get_scale_factor(units)
    x_out = x / scale_factor
    y_out = y / scale_factor

    return x_out, y_out

def utm_to_lonlat(x, y, crs):
    """Return lon, lat."""
    proj = Proj(crs)

    # scale to SI base units
    units = retrieve_units_from_crs(crs)
    scale_factor = unit_utils.get_scale_factor(units)
    x_si = x * scale_factor
    y_si = y * scale_factor

    return proj(x_si, y_si, inverse=True)

def get_utm_projection(lon, lat, units='m'):
    proj4string = get_utm_proj4string(lon, lat, units)
    return Proj(proj4string)

def get_utm_proj4string(lon, lat, units='m'):
    """Return a proj4sring."""
    zone, hem = get_utm_zone(lon, lat)
    return utm_proj4string_from_zone(zone, hem, units)

def utm_proj4string_from_zone(zone, hem, units='m'):
    """Return a proj4sring.

        Parameters
        ----------
        zone : int
        hem : 'N' or 'S'

        Returns
        -------
        proj4string : str
    """
    if hem.lower() in ('s', 'south'):
        hem_long = 'south'
    elif hem.lower() in ('n', 'north'):
        hem_long = 'north'
    else:
        raise ValueError('Unkown UTM hemisphere: %s' % hem)

    fmt = '+proj=utm +datum=WGS84 +no_defs +zone=%s +%s +units=%s'
    proj4string = fmt % (str(zone), hem_long, units)
    return proj4string


################################################################
# helpers                                                      #
################################################################
def get_utm_zone(lon, lat):
    """Return UTM zone as int an and str.

        Does not take into account any exceptions in the UTM zone definitions
        (neither those around Norway, nor those in the polar regions).

        Parameters
        ----------
        lon : float
            (deg) geographical longitude
        lat : float
            (deg) geographical latitute

        Returns
        -------
        zone : int
            a value from 1..60
        hem : 'N' or 'S'
    """
    assert -180 <= lat <= 360
    assert -90 <= lat <= 90

    # compute UTM zone
    # ------------------------------------------------
    # fraction of plant's circumference into eastward direction,
    # starting at -180W:
    fraction = ((lon + 180) / 360) % 1      # 0..1
    zone = int(fraction * 60) + 1
    # ------------------------------------------------

    # compute UTM hemisphere
    # ------------------------------------------------
    if lat > 0:
        hem = 'N'
    else:
        hem = 'S'
    # ------------------------------------------------

    return zone, hem

def retrieve_units_from_crs(s):
    """Return a str."""
    pattern = 'units='
    L = len(pattern)
    if pattern not in s:
        return 'm'

    ibeg = s.index(pattern) + L
    words = s[ibeg:].split()
    units = words[0]
    return units
