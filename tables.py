#!/usr/bin/python
"""String table utilities.

Author
------
Written in 2016
by Andreas Anhaeuser
Insitute for Geophysics and Meteorology
University of Cologne
Germany
<anhaeus@meteo.uni-koeln.de>
"""
import os
import collections

def read_vartable(filename, sep='|', comment='#', ignore_str=' '):
    """Read variable table text file and return as a dict.

        The first non-comment line is the header line. The returned dict
        consists of lists. The header line of the text file defines number and
        names of the dict keys. The subsequent lines are the entries of the
        dict, with each column grouped to its respective header.

        Parameters
        ----------
        filename : str
            location of the variable table file
        sep : str, optional
            column separator. Default: '|'
        comment : str, optional
            comment string. Set to '' or None to deactivate. Default: '#'
        ignore_str : str, optional
            The entries are shortened while they start or end on this string.
            Set to '' or None to deactivate. Default: ' '

        Returns
        -------
        table : dict
            keys are the header

        Example
        -------
        A table file that looks as shown below produces this output:
        
        {'py_name' : ['varname1', 'varname2', 'varname3'],
         'nc_name' : ['ncvn1', 'ncvn2', 'ncvn3'],
         'gain' : ['1', '1000', '1'],
         'offset' : ['0', '', '200']
         }

        py_name     | nc_name | gain | offset 
        ------------+---------+------+--------
        varname1    | ncvn1   |    1 |      0 
        varname2    | ncvn2   | 1000 |        
        varname3    | ncvn3   |    1 |    200 


        Author
        ------
        Written in 2016
        by Andreas Anhaeuser
        Institute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    if not os.path.isfile(filename):
        raise IOError('Cannot access file ' + filename)

    ###################################################
    # READ FILE                                       #
    ###################################################
    fid = open(filename, 'r')
    lines = fid.readlines()
    fid.close()

    ###################################################
    # INITIALIZE                                      #
    ###################################################
    init = False
    headers = []
    table   = {}
    S = len(sep)
    if any(ignore_str):
        I = len(ignore_str)
    else:
        I = 0
    if any(comment):
        C = len(comment)
    else:
        C = 0

    ###################################################
    # LOOP OVER LINES                                 #
    ###################################################
    for line in lines:
        ###################################################
        # SPLIT LINE INTO WORDS                           #
        ###################################################
        # the last word is removed because it is just the new line character:
        words = line.split(sep)[:-1]
        Nwords = len(words)

        # delete leading and trailing white spaces:
        for n in range(Nwords):
            word = words[n]
            while I > 0 and len(word) > 0 and word[:I]  == ignore_str:
                word = word[I:]
            while I > 0 and len(word) > 0 and word[-I:] == ignore_str:
                word = word[:-I]
            words[n] = word

        ###################################################
        # EMPTY LINE                                      #
        ###################################################
        if not any(words):
            continue

        ###################################################
        # COMMENT LINE                                    #
        ###################################################
        if C > 0 and any(words[0]) and words[0][:C] == comment:
            continue

        ###################################################
        # HEADER LINE                                     #
        ###################################################
        if not init:
            for word in words:
                table[word] = []
                headers.append(word)
            Ncols = len(headers)
            init  = True
            continue
        
        ###################################################
        # NON-DATA LINE                                   #
        ###################################################
        if len(words) != Ncols: 
            continue
        
        ###################################################
        # DATA LINE                                       #
        ###################################################
        for n in range(Ncols):
            word   = words[n]
            header = headers[n]
            table[header].append(word)

    return table

def read_namelist(filename, name_end=':', sep=',', comment='#',
        ignore_char=' ', convert_to_number=False, conversion_error=False):
    """Read namelist (setup) text file and return as a dict.

        Parameters
        ----------
        filename : str
            location of the variable table file
        name_end : str, optional
            string that marks the end of the parameter name. Default: ':'
        sep : str, optional
            string that separates elements of a value list. Default: ','
        comment : str, optional
            string that indicates the start of a comment. Default: '#'
        ignore_char : str, optional
            The entries are shortened while they start or end on this
            character(s).  Set to '' to deactivate. Default: ' '
        convert_to_number : bool, optional
            If True, entries that are not bounded by ' or " are converted to
            number (int or float). Default: False
        conversion_error : bool, optional
            If True, an error is thrown on attempted number conversion on
            invalid string. Default: False

        Returns
        -------
        table : dict
            keys are the header
        
        Example
        -------
        An ascii file like this
        
        name : Berlin
        population : 3466164
        districts : Wedding, Marzahn, Pankow        # ... and some more
        area code : 030
        # coolest clubs : outdated information
        average temperature : 9.7
        
        is returned as this dictionary:
        
        >>> namelist
        {'name' : 'Berlin',
        'population' : '3466164',
        'districts' : ['Wedding', 'Marzahn', 'Pankow'],
        'area code' : '030',
        'average temperature' : '9.7',
        }

        Author
        ------
        Written in 2016-2018
        by Andreas Anhaeuser
        Institute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    if not os.path.isfile(filename):
        raise IOError('Cannot access file ' + filename)

    ###################################################
    # READ FILE                                       #
    ###################################################
    fid = open(filename, 'r')
    lines = fid.readlines()
    fid.close()

    ###################################################
    # INITIALIZE                                      #
    ###################################################
    init = False
    namelist = {}

    E = len(name_end)
    S = len(sep)
    C = len(comment)

    ###################################################
    # LOOP OVER LINES                                 #
    ###################################################
    for line in lines:
        line = line.strip('\n')

        ###################################################
        # CUT OFF COMMENT                                 #
        ###################################################
        if C > 0 and comment in line:
            end = line.index(comment)
            line = line[:end]
        
        ###################################################
        # NON-DATA LINE                                   #
        ###################################################
        if not line.count(name_end) == 1:
            continue

        ###################################################
        # SPLIT LINE                                      #
        ###################################################
        name, data = line.split(name_end)
        name = name.strip(ignore_char)
        data = data.strip(ignore_char)

        ###################################################
        # SPLIT DATA                                      #
        ###################################################
        if S > 0 and sep in data:
            data = data.split(sep)
            for n in range(len(data)):
                data[n] = data[n].strip(ignore_char)


        ###################################################
        # FLOAT CONVERSION                                #
        ###################################################
        if convert_to_number:
            try:
                data = str2num(data)
            except ValueError as e:
                if conversion_error:
                    raise e

        ###################################################
        # WRITE TO DICT                                   #
        ###################################################
        if name in namelist.keys():
            raise KeyError(name + ' appears multiple times in namelist.')
        namelist[name] = data

    return namelist

def str2num(s):
    """Convert to number if it does not contain ' or " ."""
    # list case
    if not isinstance(s, str):
        if isinstance(s, collections.Iterable):
            return [str2num(el) for el in s]

    # check whether ' or " is in str
    str_delims = ('"', "'")
    is_string = False
    for str_delim in str_delims:
        if s[:1] == s[-1:] == str_delim:
            return s[1:-1]

    # try int
    if s.isdigit():
        number = int(s)

    # float
    else:
        number = float(s)

    return number

def column_list(headers, data, typ='s', align='l',
        column_width=16, precision=7, filename=None,
        comment_top=None, comment_bottom=None):
    """Write an ascii file with aligned columns.

        Parameters
        ----------
        headers : list of str, length N
            first line (will be commented)
        data : list of Iterable, length N
            The elements must be equal in length.
        typ : list of str, length N, optional
            {'s', 'f', 'e', 'i'} (i. e. str, float or exponential) Default: all str
        align : list of str, length N, optional
            {'l', 'c', 'r'} (i. e. left, center or right). Default: all left
        column_width : list of int, length N, optional
            Default: all 16
        precision : list of int, length N, optional
            number of digits after period. Default: all 7
        filename : str, optional
            path to file where the table is to be saved. Default: None

        Returns
        -------
        str

        Author
        ------
        Written in 2016
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    ###################################################
    # INPUT CHECK                                     #
    ###################################################

    # ========== helper functions  ======================= #
    def expand(x, J):
        """Expand variable to length J if it is singleton."""
        if isinstance(x, collections.Iterable):
            if len(x) == 1:
                return x * J
            elif len(x) == J:
                return x
            else:
                raise IndexError()
        else:
            return [x] * J

    def check_str(x):
        """Check whether variable is str or None. Raise error if not."""
        if x is None:
            return
        if isinstance(x, str):
            return
        raise TypeError()

    # ========== perform checks  ========================= #
    # headers
    assert isinstance(headers, collections.Iterable)
    J = len(headers)

    # data
    assert isinstance(data, collections.Iterable)
    assert len(data) == J
    init = False
    for x in data:
        assert isinstance(x, collections.Iterable)
        if not init:
            N = len(x)
            init = True
        assert len(x) == N

    # length J variables
    typ = expand(typ, J)
    align = expand(align, J)
    column_width = expand(column_width, J)
    precision = expand(precision, J)

    # strings
    check_str(filename)
    check_str(comment_top)
    check_str(comment_bottom)

    ###################################################
    # ABBREVIATIONS                                   #
    ###################################################
    cw = column_width
    prec = precision

    ###################################################
    # INITIALIZE                                      #
    ###################################################
    text = ''

    ###################################################
    # WRITE HEADLINES                                 #
    ###################################################
    line = ''
    for j in range(J):
        if j == 0:
            word = '# ' + headers[j]
        else:
            word = headers[j]

        # justify
        if align[j] == 'l':
            word = word.ljust(cw[j])
        elif align[j] == 'c':
            word = word.center(cw[j])
        elif align[j] == 'r':
            word = word.rjust(cw[j])

        if j > 0:
            word = ' ' + word

        # append
        line = line + word
    text = text + line.rstrip()

    ###################################################
    # FORMATTERS                                      #
    ###################################################
    fmts = []
    for j in range(J):
        # justify
        if align[j] == 'l':
            a = '<'
        elif align[j] == 'c':
            a = '^'
        elif align[j] == 'r':
            a = '>'

        # ========== digits  ============================= #
        # float
        if typ[j] in 'ef':
            # e.g. ' {:12.4f}'
            d ='%1.0f.%1.0f' % (cw[j], prec[j]) + typ[j]
        # integer
        elif typ[j] == 'i':
            raise NotImplementedError('')
            # feel free to implement this
            d = '%1.0fd' % cw[j]
        # string
        elif typ[j] == 's':
            raise NotImplementedError('')
            # feel free to implement this
            d = '%1.0fs' % cw[j]
        else:
            raise NotImplementedError('')

        # formatter
        fmt = '{:' + a + d + '}'    # e.g. ' {:<14}'
        fmts.append(fmt)

    ###################################################
    # WRITE DATA                                      #
    ###################################################
    for n in range(N):
        line = '\n'
        for j in range(J):
            word = fmts[j].format(data[j][n])
            if j > 0:
                word = ' ' + word
            line = line + word
        text = text + line.rstrip()

    ###################################################
    # INCLUDE COMMENTS                                #
    ###################################################
    if comment_top is not None:
        lines = comment_top.split('\n')
        N = len(lines)
        for n in range(N):
            lines[n] = '# ' + lines[n] + '\n'
        text = ''.join(lines) + text

    if comment_bottom is not None:
        lines = comment_bottom.split('\n')
        N = len(lines)
        for n in range(N):
            lines[n] = '\n# ' + lines[n]
        text = text + ''.join(lines)

    # add final eol:
    text = text + '\n'

    ###################################################
    # WRITE TO FILE                                   #
    ###################################################
    if filename is not None:
        idx = filename.rfind('/')
        if idx >= 0:
            path = filename[:idx]
        else:
            path = '.'
        if not os.path.isdir(path):
            os.makedirs(path)
        if os.path.isfile(filename):
            os.remove(filename)
        with open(filename, 'w') as fid:
            fid.write(text)

    return text

def get_column_list(*args, **kwargs):
    """Alias to read_column_list."""
    return read_column_list(*args, **kwargs)

def read_column_list(filename, sep=None, skip_rows=0, comment_str='#'):
    """Read text file structured in columns and return as list of lists.
    
        Parameters
        ----------
        filename : str
            path to text file
        sep : None or str, optional
            column separator, as in str.split(). Default: None
        skip_rows : int, optional
            number of rows at the beginning to skip, Detault: 0
        comment_str : str, optional
            str indicating comments. Default: '#'

        Returns
        -------
        cols : list of lists of str
            cols[i][j] contains the text of j-th row in the i-th column.
    """
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    assert isinstance(filename, str)
    assert os.path.isfile(filename)
    if sep is not None:
        assert isinstance(sep, str)
        assert len(sep) > 0
    assert isinstance(skip_rows, int)
    assert skip_rows >= 0

    ###################################################
    # READ FILE                                       #
    ###################################################
    fid = open(filename)
    lines = fid.readlines()
    fid.close()

    ###################################################
    # RETRIEVE COLUMNS                                #
    ###################################################
    init = False
    for line in lines:
        # skip_rows
        if skip_rows > 0:
            skip_rows -= 1
            continue

        # remove leading and trailing white spaces and newline:
        line = line.strip('\n').strip()
        
        # skip empty line:
        if line == '':
            continue

        # comments
        if any(comment_str):
            i = line.find(comment_str)
            # skip comment line:
            if i == 0:
                continue
            # remove in-line comment:
            if i > 0:
                line = line[:i]

        words = line.split(sep)

        # initialize output arrays
        if not init:
            cols = [[] for word in words]
            N = len(words)
            init = True

        assert len(words) == N

        for n in range(N):
            word = words[n].strip()
            cols[n].append(word)

    return cols

###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
    fn = '/home/anhaeus/py/blh/create_nc_files/comparison.table'
    cols = get_column_list(fn)
