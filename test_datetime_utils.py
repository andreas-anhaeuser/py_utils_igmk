#!/usr/bin/python3
"""Testing suite for datetime_utils."""

import datetime as dt
import unittest

import lib.pub.utils_igmk.datetime_utils as du
from lib.pub.utils_igmk.datetime_utils import Interval, DaytimePeriod, Season

class MonthsInInterval(unittest.TestCase):
    def test_months_in_interval(self):
        time_beg = dt.datetime(2000, 7, 2, 0, 3, 1, 323)
        time_end = dt.datetime(2001, 4, 5, 0, 0, 1, 11)
        time_end_flat = dt.datetime(2001, 4, 1)
        interval = Interval(time_beg, time_end)
        interval_flat = Interval(time_beg, time_end_flat)
        interval_flat_incl = Interval(
                time_beg, time_end_flat, end_inclusive=True)

        
        last_beg = dt.datetime(2001, 4, 1)
        last_end = dt.datetime(2001, 5, 1)
        last_month = Interval(last_beg, last_end)

        months = du.months_in_interval(interval)
        self.assertEqual(months[-1], last_month)

        months_flat = du.months_in_interval(interval_flat)
        self.assertNotEqual(months_flat[-1], last_month)

        months_flat_incl = du.months_in_interval(interval_flat_incl)
        self.assertEqual(months_flat_incl[-1], last_month)


if __name__ == '__main__':
    unittest.main()

