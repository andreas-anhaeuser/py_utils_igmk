#!/usr/bin/python
"""Geometry functions.

    Author
    ------
    Andreas Anhaeuser (AA)
    Insitute for Geophysics and Meteorology
    University of Cologne, Germany
    <andreas.anhaeuser@posteo.net>
"""

# standard modules
import collections
import warnings 

# PyPI modules
import numpy as np

# The spheroid values are from the World Geodetic System:
_earth_radius       = 6.3710e6    # (m) Earth's mean radius
_earth_polar_radius = 6.356752e6  # (m) Earth's polar radius
_earth_eq_radius    = 6.378137e6  # (m) Earth's equatorial radius

def distance_on_sphere(
        longitude_1,
        latitude_1,
        longitude_2,
        latitude_2,
        units='deg',
        radius=_earth_radius,
        ):
    """Return distance between two points on a speherical surface.
    
        Parameters
        ----------
        longitude_1 : float
        latitude_1 : float
        longitude_2 : float or array
        latitude_2 : float or array
        units : str, optional
            'deg' or 'rad', refers to the first four parameters, default is
            'deg'.
        sphere_radius : float or array, optional
            radius of the sphere. Default: Earth radius in metres.

        If longitude_2, latitude_2 and/or sphere_radius are given as an array,
        they must all agree in their dimensions.

        Returns
        -------
        float or array
        
        History
        -------
        2014-2016 : (AA)
    """
    ###################################
    # INPUT CHECK                     #
    ###################################
    if not isinstance(units, str):
        raise TypeError("units must be 'deg' oder 'rad'.")
    if units.lower() not in ['deg', 'degrees', 'rad', 'radians']:
        raise ValueError("units must be 'deg' oder 'rad'.")

    ###################################
    # CONVERSIONS                     #
    ###################################
    dr = np.pi / 180       # for conversion from deg to rad
    if units.lower() in ['deg', 'degrees']:
        lon1 = longitude_1 * dr
        lat1 = latitude_1  * dr
        lon2 = longitude_2 * dr
        lat2 = latitude_2  * dr
    else:
        lon1 = longitude_1
        lat1 = latitude_1
        lon2 = longitude_2
        lat2 = latitude_2

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    ###################################################
    # NOMENCLATURE                                    #
    ###################################################
    # ds : central angle (i. e. separation angle between the points)

    ########################################
    # CASE: LARGE ANGULAR SEPARATION       #
    ########################################
    # This is the algebraically exact function but it causes considerable
    # rounding errors for small angles. 
    # ds is the central angle (i. e. separation angle between the points)
    summand1 = np.sin(lat1) * np.sin(lat2)
    summand2 = np.cos(lat1) * np.cos(lat2) * np.cos(dlon)
    S = summand1 + summand2

    # numerical errors can lead to S > 1 (this is mathematically not possible,
    # however)
    if isinstance(S, collections.Iterable):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            S[S > 1] = 1.
    elif S > 1:
         S = 1.
    ds = np.arccos(S)

    ########################################
    # CASE: SMALL ANGULAR SEPARATION     #
    ########################################
    # for small central angles use haversine function (which avoids that
    # rounding errors have large influence)
    summand1 = (np.sin(dlat / 2.))**2
    summand2 = np.cos(lat1) * np.cos(lat2) * (np.sin(dlon / 2.))**2
    ds_hav = 2 * np.arcsin(np.sqrt(summand1 + summand2))
        
    ###################################################
    # REPLACE ds BY ds_small IF ds IS SMALL           #
    ###################################################
    # If ds is array-like, change appropriate elements of ds (in this case,
    # ds_is_small is a bool matrix). If ds is a scalar, then ds_is_small is a
    # bool.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        ds_is_small = ds < 1e-3
    if isinstance(ds, collections.Iterable):
        ds[ds_is_small] = ds_hav[ds_is_small]
    elif ds_is_small:
        ds = ds_hav
                
    return radius * ds

def nearest_neighbour_on_sphere(
        lon_point,
        lat_point,
        lon_array,
        lat_array,
        units='deg',
        lonlatmode='mesh',
        values=None,
        radius=_earth_radius,
        skip_nans=False,
        ):
    """Return the nearest neighbour on a sphere.
    
        lon_array and lat_array can be given in different formats, which is
        specified by lonlatmode:
        1) If you give them as a field (meshgrid), where lon[i, j] and 
           lat[i, j] are the coordinates at point [i, j], set lonlatmode to
           'mesh'. In this case, index is the pair (i, j). This works also with
           arrays that have more than 2 axes for some reason. Then index is
           returned as (i, j, k, ...)
        2) If they are a 1d-list of equal length, where lon[i] and lat[i] give
           the coordinates of point i, set lonlatmode to 'list'. In this case,
           index is the tuple (i,)
        3) If they are axes, so that the coordinates of point [i, j] are given
           by lon[i] and lat[j], then set lonlatmode to 'axes'. In this case,
           index is the pair (i, j)

        Parameters
        ----------
        lon_point : float
            longitude of the central point
        lat_point : float
            latitude of the central point
        lon_array : array of floats
            array with arbitrary number of axes
        lat_array : array of floats
            array with arbitrary number of axes
            if lonlatmode is not 'axes', then dimensions must be same as of
            lon_array
        units : 'deg' or 'rad'
            units of longitudes and latitudes
        lonlatmode : 'mesh', 'list' or 'axes'
            determines how lon_array and lat_array are interpreted.
        values : None or array or list, optional
            values at the grid points. Same length as lon_array and lat_array.
            array with arbitrary number of axes. Dimensions must be same as of
            lon_array. Default: None
        skip_nans : bool, optional
            if True, then positions in the grid are ignored if the
            corresponding element in values is nan. Default: False.
      
        Returns
        -------
        nlon :float or None
            longitude of the nearest point. None if no such point exists
        nlat : float or None
            latitude of the nearest point. None if no such point exists
        idx : tuple of ints or None
            index of nlon and nlat in lon_array and lat_array
        value : object or None
            the values[idx] (if applicable)
        dist : float
            distance of (nlon, nlat) to (lon_point, lat_point)

        History
        -------
        2014-2016 : (AA)
    """ 
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    # units:
    if not isinstance(units, str):
        raise TypeError("units must be 'deg' oder 'rad'.")
    if units.lower() not in ['deg', 'degrees', 'rad', 'radians']:
        raise ValueError("units must be 'deg' oder 'rad'.")

    # lonlatmode:
    if lonlatmode.lower() in ['m', 'mesh', 'field', '2d', 'list', 'l']:
        if not np.shape(lon_array) == np.shape(lat_array):
            raise TypeError('input arrays must be of same shape.')
        lon, lat = lon_array, lat_array
    elif lonlatmode.lower() in ['axes', 'a']:
        lon, lat = np.meshgrid(lon_array, lat_array)
    else:
        raise ValueError('unsupported value for lonlatmode.')

    ###################################################
    # TRIVIAL CASE                                    #
    ###################################################
    if not np.any(lon):
        return (None, None, None, None, np.nan)

    ###################################################
    # FIND NEAREST NEIGHBOUR                          #
    ###################################################
    D = distance_on_sphere(
            lon_point, lat_point, lon, lat, units=units,
            radius=radius)

    # exclude nan's:
    D[np.isnan(D)] = np.inf

    # exclude positions where values is nan:
    if skip_nans and np.any(values):
        D[np.isnan(values)] = np.inf

    # find nearest:
    idx = np.unravel_index(D.argmin(), D.shape)

    # check whether it is not infinity (this is the case if there is no
    # suitable point):
    if not D[idx] < np.inf:
        idx = None

    ###################################################
    # RETURN                                          #
    ###################################################
    # special case: no nearest neighbour found
    if idx is None:
        return (np.nan, np.nan, None, None, np.nan)

    # regular case
    if values is None:
        return (lon[idx], lat[idx], idx, None, D[idx])
    else:
        return (lon[idx], lat[idx], idx, values[idx], D[idx])        

def radius_on_spheroid(
        latitude, radius_e=_earth_eq_radius,
        radius_p=_earth_polar_radius):
    """Return the local radius on an oblate ellipsoid of revolution.

        Parameters
        ----------
        latitude : number (deg)
                   geographical latitude of the point.
        radius_e : number > 0
                   equatorial radius of the spheroid.
                   Default: Value for the Earth in meters.
        radius_p : number > 0
                   polar radius of the spheroid.
                   Default: Value for the Earth in meters.

        Returns
        -------
        float
            Distance of the point to the center of the spheroid.
            
        Note
        ----
        This is not equivalent to the radius of curvature at this point.

        History
        -------
        2016 : (AA)
    """

    # convert from deg to rad:
    phi = latitude * np.pi/180.

    # abbreviations:
    a = radius_e
    b = radius_p

    numerator   = (a**2 * np.cos(phi))**2 + (b**2 * np.sin(phi))**2
    denominator = (a    * np.cos(phi))**2 + (b    * np.sin(phi))**2
    return np.sqrt(numerator / denominator)

###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
#     import aa_io.aa_radiosonde as rs
#     time = dt.datetime(2011, 5, 3, 0, 0)
#     coor = rs.get_coordinate_list(time)
#     flon = np.array(coor['lon'])
#     flat = np.array(coor['lat'])
#     falt = np.array(coor['alt'])
#     lon = flon[0]
#     lat = flat[0]
    lon = lons[-1]
    lat = lats[-1]
    D = distance_on_sphere(lon, lat, lons, lats)
    flon = np.arange(-10, 10, 0.1)
    flat = np.arange(40, 50, 0.1)
    mlon, mlat = np.meshgrid(flon, flat)
    X = nearest_neighbour_on_sphere_opt(lon, lat, mlon, mlat)
    Y = nearest_neighbour_on_sphere(lon, lat, mlon, mlat)
    
   
