#!/usr/bin/python
"""A collection of functions that act on classes of the datetime module.

    Author
    ------
    Written in 2014-2018
    by Andreas Anhaeuser (AA) <andreas.anhaeuser@posteo.net>
    Institute for Geophysics and Meteorology
    University of Cologne, Germany
"""

# standard modules
import calendar
import collections
import datetime as dt

import numpy as np

#        ********************************************
#        *                  CLASSES                 *
#        ********************************************

class DaytimePeriod(object):
    """A portion of the day cycle.

        This type is unaware of its absolute (calendar) date.
        It can extend beyond midnight, e. g. 23:00 to 01:00.

        Examples
        --------
        >>> import datetime as dt
        >>> lon = 6.9
        >>> lat = 50.9
        >>> now = dt.datetime.now()
        >>> SR, SS = utc_sunrise_sunset(now, lon, lat)
        >>> day_period = DaytimePeriod(SR, SS)
        >>> if day_period.contains(now):
        >>>     print 'It is day.'
        >>> else:
        >>>     print 'It is night'

        # absolute dates do not matter:
        >>> SR_new = SR + dt.timedelta(days=100)
        >>> SS_new = SS + dt.timedelta(day=-2)
        >>> day_period = DaytimePeriod(SR_new, SS_new)
        >>> if day_period.contains(now):
        >>>     print 'It is day.'
        >>> else:
        >>>     print 'It is night'
    """

    def __init__(self, dt_start=None, dt_end=None, allow_whole_day=True):
        """dt_start and dt_end must be dt.datetime or dt.time objects.

            The absolute calendar dates of the input arguments are ignored.
            This means that calling this function with dt_start 1st Jan 1970
            23:00 and dt_end 30th Mar 1970 01:00 is equivalent to calling it
            with dt_start 1st Jan 2014 23:00 and dt_end 30th Mar 1789 01:00.
            Both will result in a time period between 23:00 and 01:00.

            Boundaries: The lower end dt_start is inclusive, the upper end
            dt_end is exclusive.

            Parameters
            ----------
            dt_start : dt.datetime or dt.time, optional
                default: midnight
            dt_end : dt.datetime or dt.time, optional
                default: midnight
            allow_whole_day : bool, optional
               applies only if start and end are equal (or equivalent) see Note
               below. Default: True.

            Note
            ----
            In case, dt_start == dt_end
            * if allow_whole_day is True, the time period will contain the
              whole day.
            * if allow_whole_day is False, the time period will not contain
              anything.

            In all other cases allow_whole_day is ignored.
        """
        #####################################
        # DEFAULT                           #
        #####################################
        if dt_start is None:
            dt_start = dt.datetime(1, 1, 1)
        if dt_end is None:
            dt_end = dt.datetime(1, 1, 1)

        #####################################
        # INPUT CHECK                       #
        #####################################
        for d in (dt_start, dt_end):
            if not (isinstance(d, dt.datetime) or isinstance(d, dt.time)):
                raise TypeError('dt_start and dt_end must be instances of ' +
                        'datetime. datetime or datetime.time.')
        if not isinstance(allow_whole_day, bool):
            raise TypeError('allow_whole_day must be a boolean.')

        ####################################
        # CONVERT TO dt.datetime           #
        ####################################
        if dt_start.__class__ is dt.time:
            start = dt.datetime.combine(dt.date(1, 1, 1), dt_start)
        else:
            start = dt_start
        if dt_end.__class__ is dt.time:
            end = dt.datetime.combine(dt.date(1, 1, 1), dt_end)
        else:
            end = dt_end

        ####################################
        # SHIFT TO YEAR 1                  #
        ####################################
        # self.start will be on Jan 1st 0001 CE.
        # self.end will be between 0 and 1 day later than self.end, all will
        # thus be on Jan 1st or Jan 2nd in year 0001 CE.
        start = start.replace(1, 1, 1)
        end = end.replace(1, 1, 1)

        ####################################
        # CHECK SEQUENCE                   #
        ####################################
        # make sure that end is not earlier that start:
        if end < start:
            end = end.replace(day=2)
        if end == start and allow_whole_day:
            end = end.replace(day=2)

        ####################################
        # CREATE FIELDS                    #
        ####################################
        self.start = start
        self.end = end

    def __repr__(self):
        return str(self)

    def __str__(self):
        s = self.start
        e = self.end
        text = ('DaytimePeriod from %02.0f:%02.0f:%02.0f to ' +
            '%02.0f:%02.0f:%02.0f') % \
                ((s.hour, s.minute, s.second) + \
                 (e.hour, e.minute, e.second))
        return text

    def length(self):
        """Return an instance of datetime.timedelta."""
        return self.end - self.start

    def contains(self, d):
        """Check whether d is within DaytimePeriod or not.

        Parameters
        ----------
        d : dt.datetime or dt.time or Iterable of such

        Returns
        -------
        bool or list of such

        """
        ###################################################
        # RECURSIVELY CALL FUNCTION                       #
        ###################################################
        if isinstance(d, collections.Iterable):
            return [self.contains(dd) for dd in d]

        ####################################
        # INPUT CHECK                      #
        ####################################
        if d.__class__ not in [dt.datetime, dt.time]:
            raise TypeError('Argument must be an instance of ' +
                    'datetime. datetime or datetime.time.')

        ####################################
        # CONVERT TO dt.datetime           #
        ####################################
        if d.__class__ is dt.time:
            dd = dt.datetime.combine(dt.date(1, 1, 1), d)
        else:
            dd = d.replace(1, 1, 1)

        ####################################
        # CHECK RELATIVE POSITIONS         #
        ####################################
        # make sure that dd is later than self.start:
        if dd < self.start:
            dd = dd.replace(day=2)
        # check whether dd is earlier than self.end:
        if dd < self.end:
            return True
        else:
            return False

class Season(object):
    """A section of the year cycle.

        This type is unaware of its absolute year number. It can extend beyond
        New Year, e. g. Nov 2nd to Jan 10th.

        Examples
        --------
        >>> import datetime as dt
        >>> lon = 6.9
        >>> lat = 50.9
        >>> now = dt.datetime.now()
        >>> beg = dt.datetime(1, 6, 21)
        >>> end = dt.datetime(1, 9, 23)

        # instantiate with beg and end:
        >>> season = Season(beg, end)
        >>> if season.contains(now):
        >>>     print 'It is summer!'
        >>> else:
        >>>     print 'It is not summer.'

        # alternatively, instantiate with months:
        >>> season = Season(months='DJF')
        >>> if season.contains(now):
        >>>     print 'It is meteorological winter.'
        >>> else:
        >>>     print 'It is not meteorological winter'

        # show length of season:
        >>> print season.length()

        # special case: beg == end
        >>> season = Season(beg, beg)
        >>> print season.length()
        365 days, 0:00:00
        >>> season = Season(beg, beg, allow_whole_year=False)
        >>> print season.length()
        0:00:00

        # the function is unaware of the absolute year, so pay attention in
        # leap years!
        >>> beg = dt.datetime(2016, 2, 1)
        >>> end = dt.datetime(2016, 3, 1)
        >>> season = Season(beg, end)
        >>> print season.length()
        28 days, 0:00:00

        # the function is microsecond-precise:
        >>> beg = dt.datetime(1, 6, 1)
        >>> end = dt.datetime(1, 6, 30, 12)
        >>> instant1 = end - dt.timedelta(microseconds=2)
        >>> instant2 = end + dt.timedelta(microseconds=2)
        >>> season = Season(beg, end)
        >>> print season.contains(instant1)
        True
        >>> print season.contains(instant2)
        False
       """
    def __init__(
            self, dt_start=dt.datetime(1, 1, 1), dt_end=dt.datetime(1, 1, 1),
            months='', allow_whole_year=True
        ):
        """dt_start and dt_end must be dt.datetime or dt.date objects.

            The absolute year numbers of the input arguments are ignored. This
            means that calling this function with dt_start 1st Jan 1970 and
            dt_end

            30th Mar 1970 is equivalent to calling it with dt_start 1st Jan
            2014 and dt_end 30th Mar 1789.

            Parameters
            ----------
            dt_start : dt.datetime, optional
                lower boundary, inclusive
            dt_end : dt.datetime, optional
                upper boundary, exclusive
            months : str, optional
                This overrides `dt_start` and `dt_end`. Valid values are:
                1. three-char abbreviation of one month, such as
                    {'jan', 'feb', ...} or
                2. sequence of chars of consecutive month initials, such as
                    {'djf', 'mam', ...}, but also {'fm', 'jjaso', ...}
                    (first and last months are inclusive) or
                3. 'year' for the whole year or
                4. '' to disable and use `dt_start` and `dt_end` instead
                    (default)
                `months` is not case-sensitive
            allow_whole_year : bool, optional
                only applies if `dt_start` and `dt_end` are identical and, in
                this case, determines whether Season contains the whole year or
                nothing.

            Special cases
            -------------
            1. dt_start == dt_end
                * if allow_whole_year : Season contains the whole year.
                * else                : Season does not contain anything.

            2. February 29th
                If the date of dt_start and/or dt_end is Feb 29th, they are
                treated as Mar 1st 00:00.
        """
        ###################################################
        # INPUT CHECK                                     #
        ###################################################
        for d in [dt_start, dt_end]:
            assert isinstance(d, dt.date)
        assert isinstance(allow_whole_year, bool)

        # months:
        if months.lower() == 'year':
            mon = 'jfmamjjasond'
        else:
            mon = months.lower()
        allowed_month_seqs = 'jfmamjjasondjfmamjjasond'
        allowed_months = [
                'jan', 'feb', 'mar', 'apr', 'may', 'jun',
                'jul', 'aug', 'sep', 'oct', 'nov', 'dec',
                ]
        if not isinstance(mon, str):
            raise TypeError('months must be a string.')
        if len(mon) == 1:
            raise ValueError('months must be a string of at least two ' +
                'month initial letters.')
        if len(mon) > 12:
            raise ValueError('months must not contain more than 12 letters.')
        if (mon not in allowed_months and
                mon[:3].lower() not in allowed_month_seqs):
            raise ValueError('months is not an allowed sequence of ' +
                'month initial letters.')

        ###################################################
        # INITIALIZE                                      #
        ###################################################
        # self.start will be in year 1
        # self.end will be between 0 and 1 year later than self.end, and will
        # thus be in year 1 or 2
        if mon == '':
            start = year_one(dt_start)
            end   = year_one(dt_end)
        elif mon in allowed_month_seqs:
            first_month = allowed_month_seqs.find(mon) + 1
            last_month = (first_month + len(mon)) % 12
            if last_month == 0:
                last_month = 12
            start = dt.datetime(1, first_month, 1)
            end   = dt.datetime(1, last_month, 1)
        elif mon[:3].lower() in allowed_months:
            first_month = allowed_months.index(mon) + 1
            last_month  = (first_month + 1) % 12
            if last_month == 0:
                last_month = 12
            start = dt.datetime(1, first_month, 1)
            end   = dt.datetime(1, last_month, 1)


        # make sure that end is not earlier than start:
        if end < start:
            end = end.replace(year=2)
        if end == start and allow_whole_year:
            end = end.replace(year=2)

        self.start = start
        self.end   = end

    def __repr__(self):
        return str(self)

    def __str__(self):
        s = self.start.replace(year=1900)
        e = self.end.replace(year=1900)

        fmt = '%d %b'
        add = ''
        if s.hour != 0 or e.hour != 0:
            add = ' %H:%M'
        if s.second != 0 or e.second != 0:
            add = ' %H:%M:%S'
        fmt = fmt + add
        text = 'Season from ' + s.strftime(fmt) + ' to ' + e.strftime(fmt)
        return text
    
    def months(self):
        """Return a string in the form of 'JJA' of 'Jun'.
            
            A month is considered to be 'covered' if at least 15 of its days
            are within the season.

            If the season covers more than one month, a sequence of capitals,
            such as 'JJA' (June-August) or 'NDJFMA' (November-April) is
            returned.

            If the season covers one month, than the first three characters of
            its name are returned.

            If the season covers less than 15 days in any month, '' is
            returned.
        """
        raise NotImplementedError('')

    def length(self):
        """Return an instance of datetime.timedelta."""
        return self.end - self.start

    def contains(self, d):
        """Check whether d is within Season or not.

            Parameters
            ----------
            d : dt.datetime or dt.date or Iterable of such

            Returns
            -------
            bool or list of such
        """
        ###################################################
        # RECURSIVELY CALL FUNCTION                       #
        ###################################################
        if isinstance(d, collections.Iterable):
            return [self.contains(dd) for dd in d]

        ###################################################
        # INPUT CHECK                                     #
        ###################################################
        if not isinstance(d, dt.date):
            raise TypeError('Argument must be datetime.datetime or ' + \
                    'datetime.date.')

        dd = year_one(d)
        if dd < self.start:
            dd = dd.replace(year=2)

        # check whether dd is earlier than self.end:
        if dd < self.end:
            return True
        else:
            return False

class Interval(object):
    """A time interval.

            Parameters
            ----------
            start : datetime.datetime
            end : datetime.datetime
            start_inclusive : bool, optional
                Default: True
            end_inclusive : bool, optional
                Default: False
        """
    def __init__(self, start, end, start_inclusive=True, end_inclusive=False):
        # input check
        if not isinstance(start, dt.datetime):
            raise TypeError('`start` must be datetime.datetime.')
        if not isinstance(end, dt.datetime):
            raise TypeError('`end` must be datetime.datetime.')

        self.start = start
        self.end = end

        self.start_inclusive = start_inclusive
        self.end_inclusive = end_inclusive

    def __str__(self):
        """Return a str."""
        if self.start_inclusive:
            lower_par = '['
        else:
            lower_par = '('
        if self.end_inclusive:
            upper_par = ']'
        else:
            upper_par = ')'

        lower_bound = str(self.start)
        upper_bound = str(self.end)

        s = '%s%s, %s%s' % (lower_par, lower_bound, upper_bound, upper_par)
        return s

    def __repr__(self):
        """Return a str."""
        return 'Interval %s' % str(self)

    def length(self):
        """Return a datetime.timedelta."""
        return self.end - self.start

    def contains(self, time):
        """Return a bool.

            Check whether `time` is in the interval.

            Parameters
            ----------
            time : datetime.datetime

            Returns
            -------
            bool
        """
        if not isinstance(time, dt.datetime):
            raise TypeError('`time` must be datetime.datetime.')

        if self.start_inclusive:
            lower_cond = self.start <= time
        else:
            lower_cond = self.start < time

        if self.end_inclusive:
            upper_cond = time <= self.end
        else:
            upper_cond = time < self.end

        return (lower_cond and upper_cond)


#        ********************************************
#        *      helper functions to classes         *
#        ********************************************

def year_one(d):
    """Return a datetime.datetime in year 1 CE.

        This is designed as a helper function for the Season class.

        Parameters
        ----------
        d : an instance of datetime.datetime or datetime.date

        Returns
        -------
        A datetime.datetime object. Date and time is be the same as of the
        input object, only the year is set to 1.

        Note
        ----
        Dates on 29th Feb will are converted to 1st Mar 00:00.

        Tested
        ------
        Intensively. Bug-free.
    """
    # special case: Feb 29th:
    if d.month == 2 and d.day == 29:
        return dt.datetime(1, 3, 1)
    # datetime.date object:
    if d.__class__ is dt.date:
        return dt.datetime.combine(d.replace(year=1), dt.time())
    # datetime.datetime object:
    if d.__class__ is dt.datetime:
        return d.replace(year=1)

def daytimediff(minuend, subtrahend, mode=None):
    """Return the difference in daytime regardless of the absolute date.

        Parameters
        ----------
        minuend, subtrahend : datetime.datime or datetime.time objects

        Returns
        -------
        datetime.timedelta

        mode: (None, 'abs', 'pos', 'neg')
        * None: a value between -12h and + 12h is returned
        * 'abs': a value between 0 and 12h is returned. The absolute difference
        * 'pos': a value between 0 and 24h is returned.
        * 'neg': a value between -24h and 0 is returned.

        Example
        -------
        >>> minuend    = datetime.datetime(2015, 1, 1, 12, 0, 0)
        >>> subtrahend = datetime.datetime(1970,12,16,  8, 0, 0)
        >>> diff = daytimediff(minuend, subtrahend)
        >>> print(repr(diff))
        datetime.timedelta(0, 14400)
        >>> print diff
        4:00:00

        Reliability
        -----------
        Moderately tested. Should be bug-free.
    """
    ###################################
    # RECURSIVELY CALL FUNCTION       #
    ###################################
    if isinstance(minuend, list):
        return [daytimediff(m, subtrahend) for m in minuend]
    if isinstance(subtrahend, list):
        return [daytimediff(minuend, s) for s in subtrahend]

    ###################################
    # INPUT CHECK                     #
    ###################################
    for d in [minuend, subtrahend]:
        if d.__class__ not in [dt.datetime, dt.time]:
            raise TypeError(
          'Arguments must be instances of datetime.datime or datetime.time.')

    if isinstance(minuend, dt.datetime):
        m = minuend.replace(1, 1, 1)
    else:
        hour = minuend.hour
        minute = minuend.minute
        second = minuend.second
        microsec = minuend.microsecond
        m = dt.datetime(1, 1, 1, hour, minute, second, microsec)

    if isinstance(subtrahend, dt.datetime):
        s = subtrahend.replace(1, 1, 1)
    else:
        hour = subtrahend.hour
        minute = subtrahend.minute
        second = subtrahend.second
        microsec = subtrahend.microsecond
        s = dt.datetime(1, 1, 1, hour, minute, second, microsec)

    diff = m - s
    if mode is None:
        while diff.total_seconds() > 86400/2:
            diff -= dt.timedelta(days=1)
        while diff.total_seconds() < -86400/2:
            diff += dt.timedelta(days=1)

    if mode == 'abs':
        diff = abs(diff, -diff)

    if mode == 'pos':
        while diff.total_seconds() < 0:
            diff += dt.timedelta(days=1)

    if mode == 'neg':
        while diff.total_seconds() > 0:
            diff -= dt.timedelta(days=1)

    return diff


#        ********************************************
#        *                 FUNCTIONS                *
#        ********************************************

###################################################
# range                                           #
###################################################
def date_range(beg, end, inc=None, season_of_year=None):
    """Return a list of dt.date. Behaves in analogy to datetime_range().

        Works as datetime_range(). Refer there for documentation. 

        Parameters
        ----------
        beg : datetime.date
            inclusive
        end : datetime.date
            exclusive
        inc : datetime.timedelta, optional
            increment. Must be a multiple of 1 day. Default: 1 day
        season_of_year: Season, optional
            Default: whole year

        Returns
        -------
        list of datetime.date

        History
        -------
        2018-12-20 (AA): Force `inc` to be a multiple of 1 day
        2018-01-04 (AA): Made `inc` an optional argument (default: 1 day).
        2017-12-30 (AA): Created
    """
    # ========== default ================================= #
    if inc is None:
        inc = dt.timedelta(days=1)

    # ========== input check ============================= #
    # make sure inc is multiple of a day
    if inc % dt.timedelta(days=1) != dt.timedelta():
        raise ValueError('`inc` must be a multiple of 1 day.')

    # ========== normalize bounds ======================== #
    # make them datetime.datetime objects

    # the start is always midnight
    time_beg = dt.datetime.combine(beg, dt.time())

    # the end is set to midnight only if it is not already a dt.datetime
    if isinstance(end, dt.datetime):
        time_end = end
    else:
        time_end = dt.datetime.combine(end, dt.time())

    # ========== construct the list ====================== #
    dtrange = datetime_range(
            beg=time_beg, end=time_end, inc=inc,
            season_of_year=season_of_year)

    # ========== re-cast ================================= #
    # convert to datetime.date objects
    return [t.date() for t in dtrange]

def datetime_range(
        beg,
        end,
        inc,
        season_of_year=None,
        daytime_period=None,
        ):
    """Return a list of dt.datetime. Works in analogy to range().

        The function works similar to the standard range() function and is
        intended to be an extension of it to datetime objects.

        The returned list contains only elements that lie within
        season_of_year and daytime_period.

        Parameters
        ----------
        beg : dt.datetime
            inclusive
        end : dt.datetime
            exclusive
        inc : dt.timedelta
        season_of_year : Season, optional
            Default: whole year
        daytime_period : DaytimePeriod, optiona
            Default : whole day

        Tested
        ------
        Extensively tested and heaviliy used. Should be bug-free.

        History
        -------
        2014       (AA): Created
        2017-12-30 (AA): Extention for inc < 0. Raise error on inc == 0.
    """
    ###################################################
    # DEFAULT                                         #
    ###################################################
    if season_of_year is None:
        season_of_year = Season(dt.datetime(1, 1, 1), dt.datetime(1, 1, 1))
    if daytime_period is None:
        daytime_period = DaytimePeriod(dt.time(), dt.time())

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    if not isinstance(beg, dt.datetime):
        raise TypeError('start must be datetime.datetime.')

    if not isinstance(end, dt.datetime):
        raise TypeError('end must be datetime.datetime.')

    if not isinstance(inc, dt.timedelta):
        raise TypeError('inc must be datetime.timedelta.')
    if inc == dt.timedelta():
        raise ValueError('inc must not be 0.')

    if not isinstance(season_of_year, Season):
        raise TypeError('season_of_year must be Season.')

    if not isinstance(daytime_period, DaytimePeriod):
        raise TypeError('daytime_period must be DaytimePeriod.')

    ###################################################
    # ABBREVIATIONS                                   #
    ###################################################
    season = season_of_year
    timeperiod = daytime_period

    ###################################################
    # BUILD LIST                                      #
    ###################################################
    # case: positive increment
    if inc > dt.timedelta():
        proceed = lambda d : d < end
    # case: negative increment
    else:
        proceed = lambda d : d > end

    out = []
    d = beg
    while proceed(d):
        if season.contains(d) and timeperiod.contains(d):
            out.append(d)
        d += inc
    return out

###################################################
# unixtime (seconds since ...)                    #
###################################################
def datetime_to_seconds(d, reference_date=dt.datetime(1970, 1, 1)):
    """Convert datetime.datetime to seconds since reference date.

        None values are converted to nan's.

        The function also works for d given as nested lists, but only if the
        length of the sub-lists are equal (i. e. if it is matrix-like).

        Parameters
        ----------
        d : datetime.datetime or datetime.date or None, or list of such
        reference_date : datetime.datetime

        Returns
        -------
        seconds : float or array of such
            (s) time since refernce date

        Development status
        ------------------
        Extensively tested and heavily used. Should be bug-free.

        Note
        ----
        This function is intended to be 100% consistent with its reverse
        function `seconds_to_datetime`. If you find exceptions to this, the
        author would be very thankful for a note!

        See also
        --------
        datetime_to_seconds : reverse function
        """
    ref = reference_date

    # convert date to datetime:
    if ref.__class__ == dt.date:
        ref = dt.datetime(ref.year, ref.month, ref.day)

    if d.__class__ == dt.date:
        d = dt.datetime(d.year, d.month, d.day)

    # recursively call function:
    if isinstance(d, collections.Iterable):
        if len(d) == 0:
            return []
        else:
            s = [datetime_to_seconds(dd, reference_date=ref) for dd in d]
            return np.array(s)

    if d is None:
        return np.nan
    else:
        distance = d - ref                 # (this is an instance of timedelta)
        return distance.total_seconds()

def seconds_to_datetime(seconds, reference_date=dt.datetime(1970, 1, 1)):
    """Convert seconds since reference date to dt.datetime.

        nan's are converted to None-values.

        Parameters
        ----------
        seconds : float or array of such
        reference_date : datetime.datetime

        Returns
        -------
        time : datetime.datetime or None, or list of such
            If the input is given in form of array, the output is a list of
            corresponding shape

        Development status
        ------------------
        Extensively tested and heavily used. Should be bug-free.

        Note
        ----
        This function is intended to be 100% consistent with its reverse
        function `datetime_to_seconds`. If you find exceptions to this, the
        author would be very thankful for a note!

        See also
        --------
        datetime_to_seconds : reverse function
    """
    ref = reference_date

    # convert date to datetime:
    if ref.__class__ == dt.date:
        ref = dt.datetime(ref.year, ref.month, ref.day)

    if isinstance(seconds, collections.Iterable):
        return [seconds_to_datetime(s, ref) for s in seconds]

    if np.isnan(seconds):
        return None
    else:
        return ref + dt.timedelta(seconds=float(seconds))

###################################################
# JULIAN DAYS                                     #
###################################################
def julian_days_to_datetime(days):
    """Return a dt.datetime.

        Parameters
        ----------
        days : float or an Iterable of such

        Returns
        -------
        dt.datetime or list of such

        Tested
        ------
        Moderately. Seems to be working but no guarantee.
    """
    if isinstance(days, collections.Iterable):
        return([julian_days_to_datetime(d) for d in days])
    JD0 = 1721425.5  # Julian date of 1st Jan 1, 00:00
    return dt.datetime(1, 1, 1) + dt.timedelta(days=days - JD0)

def datetime_to_julian_days(time):
    """Return a float.

        Parameters
        ----------
        time : datetime.datetime object or list of such

        Returns
        -------
        float or array of such
    """
    if isinstance(time, collections.Iterable):
        return np.array([datetime_to_julian_days(t) for t in time])
    JD0 = 1721425.5  # Julian date of 1st Jan 1, 00:00
    diff = (time - dt.datetime(1, 1, 1))
    D = diff.days
    S = diff.seconds / 86400.
    U = diff.microseconds / (86400 * 1e6)
    return D + S + U + JD0

###################################################
# DAY OF YEAR                                     #
###################################################
def doy(date):
    """Return day of year and year as pair of ints.

        Deprecated. Use sod instead!

        Parameters
        ----------
        date : datetime.date or datetime.datetime or list of such

        Returns
        -------
        doy : int or array of such
            day of year. Between 1 and 366 (inclusive, in leap years)
        year : int or array of such

        Note
        ----
        1st January is DOY 1 (not 0).
    """
    ###################################################
    # RECURSIVE FUNCTION CALL FOR LISTS               #
    ###################################################
    if isinstance(date, collections.Iterable):
        S = np.shape(date)
        assert len(S) == 1
        N = S[0]

        days = np.zeros(N, dtype=int)
        years = np.zeros(N, dtype=int)
        for n in range(N):
            days[n], years[n] = doy(date[n])
        return days, years

    ###################################################
    # TYPE CONVERSION                                 #
    ###################################################
    if isinstance(date, dt.datetime):
        date = date.date()

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    assert isinstance(date, dt.date)

    year = date.year
    ref = dt.date(year, 1, 1)
    day = (date - ref).days + 1
    return day, year

def date_from_doy(doy, year):
    """Return a datetime.date.

        Parameters
        ----------
        doy : int
            day of year. Can be negative or positve. Can be larger than 366.
        year : int

        Returns
        -------
        date : datetime.date or datetime.datetime

        Note
        ----
        1st January is DOY 1 (not 0).
    """
    return dt.date(year, 1, 1) + dt.timedelta(days=doy-1)

def datetime_from_doy(doy, year):
    """Return a datetime.date.

        Parameters
        ----------
        doy : int
            day of year. Can be negative or positve. Can be larger than 366.
        year : int

        Returns
        -------
        date : datetime.date or datetime.datetime

        Note
        ----
        1st January is DOY 1 (not 0).
    """
    return dt.datetime(year, 1, 1) + dt.timedelta(days=doy-1)

def days_in_month(date):
    """Return the number of days of that month.
    
        Parameters
        ----------
        date : dt.date
        
        Returns
        -------
        int
            number of days of that months.
    """
    year = date.year
    month = date.month
    if month in [1, 3, 5, 7, 8, 10, 12]:
        return 31
    elif month in [4, 6, 9, 11]:
        return 30
    elif calendar.isleap(year):
        return 29
    else:
        return 28

###################################################
# SECONDS OF DAY                                  #
###################################################
def sod_doy(time):
    """Return second of day, day of year and year as triple of ints.

        Parameters
        ----------
        time : datetime.datetime or list of such

        Returns
        -------
        sod : int or array of such
            second of day. Between 1 and 86400 (inclusive)
        doy : int or array of such
            day of year. Between 1 and 366 (inclusive, in leap years)
        year : int or array of such

        Notes
        -----
        1st January is DOY 1 (not 0).
        00:00:00 is SOD 1 (not 0).
    """
    ###################################################
    # RECURSIVE FUNCTION CALL FOR LISTS               #
    ###################################################
    if isinstance(time, collections.Iterable):
        S = np.shape(time)
        assert len(S) == 1
        N = S[0]

        secs = np.zeros(N, dtype=int)
        days = np.zeros(N, dtype=int)
        years = np.zeros(N, dtype=int)
        for n in range(N):
            secs[n], days[n], years[n] = sod_doy(time[n])
        return secs, days, years

    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    assert isinstance(time, dt.datetime)

    year = time.year
    month = time.month
    day = time.day

    ref_doy = dt.datetime(year, 1, 1)
    ref_sod = dt.datetime(year, month, day)
    day_of_year = (time - ref_doy).days + 1
    sec_of_day = (time - ref_sod).seconds + 1
    return sec_of_day, day_of_year, year

###################################################
# ITERATION                                       #
###################################################
def next_month(year, month):
    """Return (year, month) of the following month.
    
        Parameters
        ----------
        year : int
        month : int
        
        Returns
        -------
        year_next : int
        month_next : int
    """
    if month < 12:
        month += 1
    else:
        month = 1
        year += 1
    return year, month

def next_day(year, month, day):
    """Return (year, month, day) of the following day.
    
        Parameters
        ----------
        year : int
        month : int
        day : int
        
        Returns
        -------
        year_next : int
        month_next : int
        day_next : int
    """
    thisday = dt.datetime(year, month, day)
    nextday = thisday + dt.timedelta(days=1)
    y = nextday.year
    m = nextday.month
    d = nextday.day
    return y, m, d

###################################################
# STRINGS                                         #
###################################################
def name_of_month(month, kind='full'):
    """Return a str which represents the name of the month.

        Parameters
        ----------
        month : int, between 1 and 12.
            number of the month
        kind : {'full', 'short', 'abbr'}, optional
            'full' : return full month name
            'short' : return the first 3-characters of the month name with no
                      trailing '.'
            'abbr' : return a max. 4-char str, Either the month name, if it is
                     no longer than 4, otherwise its first 3 characters plus
                     trailing '.'
            Default : 'full'

        Returns
        -------
        str
            The month name.
    """
    assert isinstance(month, int)
    assert 0 < month < 13

    names = ['January', 'February', 'March', 'April', 'May', 'June',
            'July', 'August', 'September', 'October', 'November',
            'December']

    if kind == 'full':
        pass
    elif kind == 'short':
        # shorten names
        names = [name[:3] for name in names]
    elif kind == 'abbr':
        # abbreviate names
        for n, name in enumerate(names):
            if len(name) > 4:
                names[n] = name[:3] + '.'
    else:
        raise ValueError('Unrecognized kind: ' + str(kind))

    return names[month-1]


#        ********************************************
#        *                  TESTING                 *
#        ********************************************
if __name__ == "__main__":
    beg = dt.datetime(2017, 1, 1)
    end = dt.datetime(2017, 1, 2, 12, 3)

    t1 = dt.datetime(2017, 1, 2, 2)
    t2 = dt.datetime(2017, 1, 4, 2)
    t3 = dt.datetime(2016, 1, 2, 2)

    i = Interval(beg, end)
    print(i)
    print(i.length())
    print(i.contains(t1))
    print(i.contains(t2))
    print(i.contains(t3))
    print(repr(i))

