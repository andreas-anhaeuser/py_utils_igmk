#!/usr/bin/python3
"""Testing suite for chronometer.Chronometer."""

import unittest
from time import sleep

from lib.pub.utils_igmk.chronometer import Chronometer

###################################################
# TESTING                                         #
###################################################
class BasicFunctionality(unittest.TestCase):
    def setUp(self):
        pass

    def test_loop_and_show(self):
        N = 2000
        chrono = Chronometer(N)
        for i in range(N):
            if i % 100 == 0:
                chrono.issue(i)
            sleep(0.005)
            chrono.loop_and_show()

        chrono.resumee()

    def test_init_with_iterable(self):
        items = '*' * 2000
        header = 'Initialized with Iterable'
        chrono = Chronometer(items, header=header)
        N = len(items)
        for i in range(N):
            if i % 100 == 0:
                chrono.issue(i)
            sleep(0.005)
            chrono.loop_and_show()

        chrono.resumee()


if __name__ == '__main__':
    unittest.main()
