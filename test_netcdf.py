#!/usr/bin/python3
"""Testing suite for datetime_utils."""

import unittest

import utils_igmk.netcdf as nc

class SelectVariables(unittest.TestCase):
    def test_select(self):
        filename_in = 'test_files/file_in.nc'
        filename_out = 'test_files/file_out.nc'
        varnames = ('time', 'qa_value')
        nc.select_variables(filename_in, filename_out, varnames)

        print(
                'Selected %s from\n%s -> %s.\nCheck manually.'
                % (varnames, filename_in, filename_out)
                )

if __name__ == '__main__':
    unittest.main()
