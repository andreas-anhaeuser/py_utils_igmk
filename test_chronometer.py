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
        N = 1000
        header = 'Initialized with integer'
        chrono = Chronometer(N, header=header)
        for i in range(N):
            if i % 100 == 0:
                chrono.issue(i)
            sleep(0.005)
            chrono.loop_and_show()

        chrono.resumee()

    def test_init_with_iterable(self):
        N = 1000
        items = '*' * N
        header = 'Initialized with Iterable'
        chrono = Chronometer(items, header=header)
        self.assertEqual(chrono.get_total_count(), N)
        for i in range(N):
            if i % 100 == 0:
                chrono.issue(i)
            sleep(0.005)
            chrono.loop_and_show()

        chrono.resumee()

class Colors(unittest.TestCase):
    def setUp(self):
        self.N = 300
        self.sleep_time = 0.005
        self.chrono = Chronometer(self.N)

    def sleep(self):
        sleep(self.sleep_time)

    def test_color(self):
        chrono = self.chrono
        chrono.header = 'With colors'
        for i in range(self.N):
            if i % 100 == 0:
                chrono.issue(i)
            chrono.loop_and_show()
            self.sleep()
        chrono.resumee()

    def test_no_color(self):
        chrono = self.chrono
        chrono.header = 'Without colors'
        chrono.print_colors = False
        for i in range(self.N):
            if i % 100 == 0:
                chrono.issue(i)
            chrono.loop_and_show()
            self.sleep()
        chrono.resumee()

if __name__ == '__main__':
    unittest.main()
