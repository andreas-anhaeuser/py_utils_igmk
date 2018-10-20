#!/usr/bin/python
"""Functions related to Sun's position on the sky.

    Author
    ------
    2014-2018
    Andreas Anhaeuser (AA) <andreas.anhaeuser@posteo.net>
    Institute for Geophysics and Meteorology
    University of Cologne, Germany
"""
import ephem

def utc_sunrise_sunset(d, lon, lat, alt=0, pres=None, temp=None):
    """Return SR and SS for this day as a pair of dt.datetime.

        The time of day of the input is ignored, only the date is taken into
        account.

        The horizon is assumed to be at 0 degrees, regardless of the altitude.

        Parameters
        ----------
        d : datetime.date or datetime.datetime
        lon: float
            the observer's longitude in deg North
        lat: float
            the observer's latitude in deg East
        alt: float
            the observer's altitude above the surface in meters
        pres : float, optional
            (Pa) pressure at ground. Used to compute atmospheric light
            refraction. Note: overrides alt !
        temp : float, optional
            (K) temperature at ground. Used to compute atmospheric light
            refraction. Default: 288.15.

        Returns
        -------
        SR, SS : datetime.datetime or None
            SR and SS denote sun rise and sun set of the day specified by d.
            If (in polar regions,) there is no sun set or sun rise at this day,
            SR and/or SS are None.

        Notes
        -----
        All times and dates are treated as UTC (not local time). This means
        that SS can be before SR (this is roughly the case for places that lie
        on the hemisphere where abs(lon) > 90).

        SR and SS as returned by this function are always on the same UTC date
        as d. If such a SR and/or SS does not exist, then the respective value
        is set to None. This may be counter-intuitive in some cases, e. g. on a
        day where the previous sun rise is just before 0:00 and the following
        is just after 0:00 (i. e. 24 h and some minutes later), then on the
        actual day, there is no sun rise and None is returned:

              | day (n-1) | day n   | day (n+1)
        ------+-----------+---------+-----------
        SR    | 23:59     | ----    | 00:01
        SS    | 11:20     | 11:18   | 11:16

        Dependencies
        ------------
        Makes use of the package ephem.

        Raises
        ------
        TypeError, ValueError

        Tested
        ------
        Moderately tested. Seems bug-free, but not 100% sure.
    """
    #############################################
    # INPUT CHECK                               #
    #############################################
    if not isinstance(d, dt.date):
        raise TypeError('d must be an instance of datetime.date.')
    if not -180 <= lon <= 180:
        raise ValueError('lon must be between -180 and 180.')
    if not -90 <= lat <= 90:
        raise ValueError('lat must be between -90 and 90.')

    #############################################
    # SET TIME TO MIDNIGHT                      #
    #############################################
    if isinstance(d, dt.datetime):
        d = d.date()

    #############################################
    # CREATE EPHEM OBJECTS                      #
    #############################################
    # create Sun object:
    sun = ephem.Sun(d)

    # create observer object:
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2  # (convert from Pa to hPa)
    if temp is not None:
        site.temp = temp - 273.15    # (convert from deg C to K)
    site.date = d

    #############################################
    # SR AND SS                                 #
    #############################################
    try:
        SR = site.next_rising(sun).datetime()
        # make sure SR is on the same day:
        if SR.date() != d:
            SR = None
    except ephem.NeverUpError:
        SR = None
    except ephem.AlwaysUpError:
        SR = None

    try:
        SS = site.next_setting(sun).datetime()
        # make sure SS is on the same day:
        if SS.date() != d:
            SS = None
    except ephem.NeverUpError:
        SS = None
    except ephem.AlwaysUpError:
        SS = None

    return (SR, SS)

def is_day(d, lon, lat, alt=0., pres=None, temp=None):
    """Return a bool.

        Consistent with utc_sunrise_sunset within tens of microseconds.

        Parameters
        ----------
        d : datetime.date or datetime.datetime
            UTC
        lon: float
            the observer's longitude in deg North
        lat: float
            the observer's latitude in deg East
        alt: float
            the observer's altitude above the surface in meters
        pres : float, optional
            (Pa) pressure at ground. Used to compute atmospheric light
            refraction. Note: overrides alt !
        temp : float, optional
            (K) temperature at ground. Used to compute atmospheric light
            refraction. Default: 288.15.

        Returns
        -------
        bool

        Dependencies
        ------------
        Makes use of the package ephem.

        Raises
        ------
        TypeError, ValueError

        Tested
        ------
        Moderately tested. Seems bug-free, but not 100% sure.
    """
    #############################################
    # INPUT CHECK                               #
    #############################################
    if not isinstance(d, dt.datetime):
        raise TypeError('d must be datetime.datetime.')
    if not -180 <= lon <= 180:
        raise ValueError('lon must be between -180 and 180.')
    if not -90 <= lat <= 90:
        raise ValueError('lat must be between -90 and 90.')

    #############################################
    # CREATE EPHEM OBJECTS                      #
    #############################################
    # create Sun object:
    sun = ephem.Sun(d)

    # create observer object:
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    site.date = d
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2  # (convert from Pa to hPa)
    if temp is not None:
        site.temp = temp - 273.15    # (convert from deg C to K)

    # compute sun elevation
    sun.compute(site)
    elevation = sun.alt

    # take into account extent of the sun:
    size_arcsec = sun.size
    size = size_arcsec *np.pi / (3600 *180)
    elev_top = elevation + size / 2
    return elev_top >= 0

def is_night(d, lon, lat):
    """Return a bool.

    Description, see `is_day`."""
    return not is_day(d, lon, lat)

def last_sunrise(d, lon, lat, alt=0, pres=None, temp=None):
    """Return a dt.datetime.

        The horizon is assumed to be at 0 degrees, regardless of the altitude.

        Parameters
        ----------
        d : datetime.date or datetime.datetime
        lon: float
            the observer's longitude in deg North
        lat: float
            the observer's latitude in deg East
        alt: float
            the observer's altitude above the surface in meters
        pres : float, optional
            (Pa) pressure at ground. Used to compute atmospheric light
            refraction. Note: overrides alt !
        temp : float, optional
            (K) temperature at ground. Used to compute atmospheric light
            refraction. Default: 288.15.

        Returns
        -------
        datetime.datetime

        Notes
        -----
        Unlike utc_sunrise_sunset, this function always returns a
        datetime.datetime object. If necessary, it goes back several days until
        it finds a sunrise.

        Dependencies
        ------------
        Makes use of the package ephem.

        Tested
        ------
        Moderately tested. Seems bug-free, but not 100% sure.
    """
    # Idea
    # ====
    # 1. If sunrise on this day is earlier than d, return it.
    # 2. If sunrise is later than d, go back one day.
    # 3. In polar regions, it may be necessary to go back several days to find
    #    a sun rise, hence the while-loop.

    ###################################################
    # CREATE SUN OBJECT                               #
    ###################################################
    sun = ephem.Sun(d)

    ###################################################
    # CREATE OBSERVER OBJECT                          #
    ###################################################
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2  # Pa --> hPa
    if temp is not None:
        site.temp = temp - 273.15    # K --> deg C
    site.date = d

    ###################################################
    # FIND SUN RISE                                   #
    ###################################################
    # in extreme cases (close to the pole), the last sun rise may be up to half
    # a year ago (<= 183 days).
    d_try = d
    d_inc = dt.timedelta(days=1)
    count = 0
    found = False
    while not found:
        try:
            SR = site.previous_rising(sun).datetime()
            found = True
        except ephem.NeverUpError:
            pass
        except ephem.AlwaysUpError:
            pass
        count += 1
        d_try = d_try - d_inc
        assert count < 184        # (d) half a year
    return SR

def last_sunset(d, lon, lat, alt=0, pres=None, temp=None):
    """Return a dt.datetime.

    For description, see last_sunrise
    """
    # for comments, see last_sunrise.

    ###################################################
    # CREATE SUN OBJECT                               #
    ###################################################
    sun = ephem.Sun(d)

    ###################################################
    # CREATE OBSERVER OBJECT                          #
    ###################################################
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2
    if temp is not None:
        site.temp = temp - 273.15
    site.date = d

    ###################################################
    # FIND SUN RISE                                   #
    ###################################################
    d_try = d
    d_inc = dt.timedelta(days=1)
    count = 0
    found = False
    while not found:
        try:
            event = site.previous_setting(sun).datetime()
            found = True
        except ephem.NeverUpError:
            pass
        except ephem.AlwaysUpError:
            pass
        count += 1
        d_try = d_try - d_inc
        assert count < 184
    return event

def next_sunrise(d, lon, lat, alt=0, pres=None, temp=None):
    """Return a dt.datetime.

    For description, see last_sunrise
    """
    # for comments, see last_sunrise.

    ###################################################
    # CREATE SUN OBJECT                               #
    ###################################################
    sun = ephem.Sun(d)

    ###################################################
    # CREATE OBSERVER OBJECT                          #
    ###################################################
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2
    if temp is not None:
        site.temp = temp - 273.15
    site.date = d

    ###################################################
    # FIND SUN RISE                                   #
    ###################################################
    d_try = d
    d_inc = dt.timedelta(days=1)
    count = 0
    found = False
    while not found:
        try:
            event = site.next_rising(sun).datetime()
            found = True
        except ephem.NeverUpError:
            pass
        except ephem.AlwaysUpError:
            pass
        count += 1
        d_try = d_try - d_inc
        assert count < 184
    return event

def next_sunset(d, lon, lat, alt=0, pres=None, temp=None):
    """Return a dt.datetime.

    For description, see last_sunrise
    """
    # For comments, see last_sunrise

    ###################################################
    # CREATE SUN OBJECT                               #
    ###################################################
    sun = ephem.Sun(d)

    ###################################################
    # CREATE OBSERVER OBJECT                          #
    ###################################################
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.elevation = alt
    if pres is None:
        site.compute_pressure()
    else:
        site.pressure = pres * 1e-2
    if temp is not None:
        site.temp = temp - 273.15
    site.date = d

    ###################################################
    # FIND SUN RISE                                   #
    ###################################################
    d_try = d
    d_inc = dt.timedelta(days=1)
    count = 0
    found = False
    while not found:
        try:
            event = site.next_setting(sun).datetime()
            found = True
        except ephem.NeverUpError:
            pass
        except ephem.AlwaysUpError:
            pass
        count += 1
        d_try = d_try - d_inc
        assert count < 184
    return event

def fraction_of_day(time, sunrise, sunset):
    """Return a float between 0 and 1.

        Divides the day into daytime and nighttime. Daytime ranges from 0 to
        0.5, nighttime from 0.5 to 1.0. If day and night do not have the same
        (physical) length, then the mapping from time to day_fraction is NOT
        linear, i. e. the interval from day_fraction 0.4 to 0.5 is longer or
        shorter (in seconds) than the interval from 0.5 to 0.6.

        Parameters
        ----------
        time : datetime.datetime
        sunrise : datetime.datetime
        sunset : datetime.datetime

        Returns
        -------
        float
            0. : sun rise
            0.5 : sun set
            1. sun rise
    """
    for x in [time, sunrise, sunset]:
        if not isinstance(x, dt.datetime):
            raise TypeError('time, sunrise and sunset must be ' +
                    'datetime.datetime objects.')
    if not time.date() == sunrise.date() == sunset.date():
        raise ValueError('time, sunrise and sunset must be on the same day.')

    # abbreviations:
    t = time
    SR = sunrise
    SS = sunset

    #####################################
    # CASE 1: sunrise before sunset     #
    #####################################
    if SR <= SS:
        ld = (SS - SR).total_seconds()
        ln = 86400. - ld
        if t < SR:
            diff = ln - (SR - t).total_seconds()
            df = 0.5 + 0.5 * diff / ln
        elif SR <= t < SS:
            diff = (t - SR).total_seconds()
            df = 0.5 * diff / ld
        else:
            diff = (t - SS).total_seconds()
            df = 0.5 + 0.5 * diff / ln

    #####################################
    # CASE 2: sunrise after sunset      #
    #####################################
    else:
        ln = (SR - SS).total_seconds()
        ld = 86400. - ln
        if t < SS:
            diff = ld - (SS - t).total_seconds()
            df = 0.5 * diff / ln
        elif SS <= t < SR:
            diff = (t - SS).total_seconds()
            df = 0.5 + 0.5 * diff / ld
        else:
            diff = (t - SR).total_seconds()
            df = 0.5 * diff / ln
    return df

def sun_position(d, lon, lat):
    """Calculate elevation (= altitude) and azimuth.
    
        For those unfamiliar with azimuth and altitude: They describe position
        in the sky by measuring angle around the horizon, then angle above the
        horizon.
        
        Parameters
        ----------
        d : datetime.date or datetime.datetime
        lon : float
            (deg) the observer's longitude
        lat : float
            (deg) the observer's latitude
        
        Returns
        -------
        altitude : float
            (deg)
        azimuth : float
            (deg)

        Authors
        -------
        Institute for Geophysics and Meteorology
        University of Cologne, Germany
        http://www.geomet.uni-koeln.de/en/home/

        2016: Written by
        Christopher Frank <cfrank@meteo.uni-koeln.de>

        2016-09-09: Slightly modified by
        Andreas Anhaeuser <andreas.anhaeuser@posteo.net>
    """
    ###################################################
    # CREATE SUN OBJECT                               #
    ###################################################
    sun = ephem.Sun(d)

    ###################################################
    # CREATE OBSERVER OBJECT                          #
    ###################################################
    site = ephem.Observer()
    site.lon = str(lon)
    site.lat = str(lat)
    site.date = d

    sun.compute(site)
    altitude_rad = sun.alt  # radiant
    azimuth_rad = sun.az    # radiant
    altitude_deg = 180. / np.pi * altitude_rad
    azimuth_deg = 180. / np.pi * azimuth_rad
    return altitude_grad, azimuth_grad
