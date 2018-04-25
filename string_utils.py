#!/usr/bin/python
"""Functions for strings.

    Author
    ------
    Andreas Anhaeuser (AA) <andreas.anhaeuser@posteo.net>
    Institute for Geophysics and Meteorology
    University of Cologne, Germany
"""

def human_format(num, digits=2, sep='', mode='prefix', type='float'):
    """Convert large and small numbers to human with SI prefixes.
    
        Parameters
        ----------
        num : float
            the number to convert
        digits : int, optional
            total (minimum) number of decimal digits to include. Default: 2
        sep : str, optional
            separator between the number and the prefix. Default: ''
        mode : {'prefix', 'power'}, optional
            'prefix' : SI prefix
            'power' : power of 10, in tex format [ '$\\times 10^{power}$' ]
            Default: 'prefix'
        type : {'float', 'int'}, optional
            if 'int', small numbers don't have decimal digits (i. e. '4' rather
            than '4.0')

        Returns
        -------
        str

        Example
        -------
        >>> human_format(234567, 1, ' ')
        '234.6 k'
    """
    ###################################################
    # INPUT CHECK                                     #
    ###################################################
    valid_modes = ['prefix', 'power']
    valid_types = ['float', 'int']

    assert isinstance(digits, int)
    assert isinstance(sep, str)
    assert mode in valid_modes
    assert type in valid_types

    ###################################################
    # SETUP                                           #
    ###################################################
    if mode == 'prefix':
        P = ['', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
        p = ['', 'm', u'\u03bc', 'n', 'p', 'f', 'a', 'z', 'y']
    elif mode == 'power':
        powers = range(0, 30, 3)
        P = ['\\times\, 10^{%1.0f}' % power for power in powers]
        p = ['\\times\, 10^{-%1.0f}' % power for power in powers]

    ###################################################
    # MAGNITUDE                                       #
    ###################################################
    mag = 0

    # large numbers
    while abs(num) >= 1000:
        if mag >= len(P) - 1:
            break
        mag += 1
        num /= 1000.

    # small numbers
    while abs(num) < 1:
        if abs(mag) >= len(p) - 1:
            break
        mag -= 1
        num *= 1000

    ###################################################
    # NUMBER OF DIGITS                                #
    ###################################################
    # before '.'
    if abs(num) >= 100:
        n_before = 3
    elif abs(num) >= 10:
        n_before = 2
    else:
        n_before = 1

    # after '.'
    n_after = max(0, digits - n_before)

    ###################################################
    # FORMATTER                                       #
    ###################################################
    if mag == 0 and type == 'int':
        # special case: int close to 1:
        fmt = '%1.0f' + sep + '%s'
    else:
        # regular case:
        fmt = '%.' + ('%1.0f' % n_after) + 'f' + sep + '%s' 

    ###################################################
    # MATH MODE                                       #
    ###################################################
    if mode == 'power':
        fmt = '$' + fmt + '$'

    ###################################################
    # BUILD STRING                                    #
    ###################################################
    # select prefix
    if mag >= 0:
        prefix = P[mag]
    else:
        prefix = p[-mag]

    return fmt % (num, prefix)
