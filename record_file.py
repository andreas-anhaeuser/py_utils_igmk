#!/usr/bin/python3
"""Use a simple text file as a record.

    The module provide very basic functionality:
    - add a record to the file
    - check whether a record exists in the file
    - remove the record from the file
"""

# standard modules
import os

###################################################
# USER FUNCTIONS                                  #
###################################################
def add(entry, filename, allow_multiple_appearances=False):
    """Append a line to the record file.

        Creates the record file autonomously if it does not exist.

        Parameters
        ----------
        entry : str
            record to be added
        filename : str
            record file
        allow_multiple_appearances : bool, optional
            (default: False) add entry another time if it is already in the
            file.

        Returns
        -------
        added : bool
            True if `entry` was added, False otherwise
    """
    # create if not already existing
    initialize(filename, overwrite=False)

    if not allow_multiple_appearances:
        if contains(entry, filename):
            return False

    # normalize line
    line = str(entry).strip()

    # add it to file
    with open(filename, 'a') as fid:
        print(line, file=fid)

    return True

def contains(entry, filename):
    """Return a bool.

        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file

        Returns
        -------
        found : bool
            True if `entry` is in file, False otherwise
    """
    nline = find_first_appearance(entry, filename)
    if nline is None:
        return False
    else:
        return True

def remove(entry, filename, remove_multiple_appearances=True):
    """Remove record from file.

        No exception is thrown if the entry does not exists.

        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file
        remove_multiple_appearances : bool, optional
            (default: True) remove all appearances of entry if it exists
            multiple times

        Returns
        -------
        found : bool
            True if `entry` was found, False otherwise.
    """
    found = False
    repeat = True

    while repeat:
        found_once = remove_once(entry, filename)

        if found_once:
            found = True

        if not remove_multiple_appearances:
            repeat = False
        elif not contains(entry, filename):
            repeat = False

    return found

###################################################
# HELPERS                                         #
###################################################
def remove_once(entry, filename):
    """Remove first appearance of record from file.

        No exception is thrown if the entry does not exists.

        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file

        Returns
        -------
        found : bool
            True if `entry` was found, False otherwise.
    """
    nline = find_first_appearance(entry, filename)
    if nline is None:
        return False

    # read file
    with open(filename, 'r') as fid:
        lines = fid.readlines()

    # remove one line
    reduced_lines = lines[:nline] + lines[nline+1]

    # re-write file
    with open(filename, 'w') as fid:
        for line in reduced_lines:
            print(line, file=fid)

    return True

def find_first_appearance(entry, filename):
    """Return line number of first appearance.
        
        Python counting (0 == first line).

        Parameters
        ----------
        entry : str
            search item
        filename : str
            record file

        Returns
        -------
        nline : int or None
            int: line number if the entry was found
            None: file does not exist or entry not found
    """
    # file does not exist --> False
    if not os.path.isfile(filename):
        return None

    # read file
    with open(filename, 'r') as fid:
        lines = fid.readlines()

    # search contents
    hit = False
    entry = str(entry).strip()
    for nline, line in enumerate(lines):
        line = line.strip()
        if entry == line:
            hit = True
            break

    if hit:
        return nline
    else:
        return None

def initialize(filename, overwrite=True):
    """Create record file.

        Parameters
        ----------
        filename : str
            record file
        overwrite : bool, optional
            (default: True) overwrite if file alreay exists

        Returns
        -------
        None
    """
    # remove already existing file
    if os.path.isfile(filename) and overwrite:
        os.remove(filename)

    # create directory
    directory = os.path.realpath(os.path.dirname(filename))
    if not os.path.isdir(directory):
        os.makedirs(directory)

    # create file
    if not os.path.isfile(filename):
        with open(filename, 'w'):
            pass

    return None


###################################################
# TESTS                                           #
###################################################
def test_write_and_look_up():
    import datetime as dt
    filename = './some/strange/path/record_file.tmp.txt'
    entries = ['first', 2, dt.datetime.now(), ('a', 'tuple')]

    print('record file: %s' % filename)
    for entry in entries:
        print('*' * 30)

        # look up
        found = contains(entry, filename)
        print('Entry %s in record file: %s' % (entry, found))

        # add
        print('Add it.')
        add(entry, filename)

        # look up again
        found = contains(entry, filename)
        print('Entry %s in record file: %s' % (entry, found))

    print('*' * 30)
    for entry in entries:
        # look up
        found = contains(entry, filename)
        print('Entry %s in record file: %s' % (entry, found))


if __name__ == '__main__':
    test_write_and_look_up()



