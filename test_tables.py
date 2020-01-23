#!/usr/bin/python3
"""Testing suite for tables."""

# standard modules
import unittest

# PyPI modules
import numpy as np

# to be tested
import utils_igmk.tables as tab

class Table(unittest.TestCase):
    def setUp(self):
        pass

    def test_table(self):
        filename = 'tables__test_file_column_list.txt'
        sep = '|'
        convert_to_number = True

        data = tab.read_column_list_with_headers(
                filename, sep=sep, convert_to_number=convert_to_number,
                missing_str='-',
                )

        for key in sorted(data):
            print('%s: %s' % (key.ljust(32), str(data[key])))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
