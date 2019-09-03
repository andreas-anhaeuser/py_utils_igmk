"""Functions that work similar to `range` and `numpy.range`."""

import numpy as np

def logarithmic_range(
        vmin, vmax, N, add_zero=True, symmetric_about_zero=False
        ):
    """Return a logarithmic range.

        Parameters
        ----------
        vmin : float > 0
            lower bound, inclusive
        vmax : float > `vmin`
            upper bound, inclusive 
        N : int > 1
            number of element in the output range
        add_zero : bool, optional
            (default: True) make zero be the lowest level
        symmetric_about_zero : bool, optinal
            (default: False) extent to both negative and positive, with a zero
            element in the middle

        Returns
        -------
        ndarray
    """
    # cast ================================================
    vmin = float(vmin)
    vmax = float(vmax)
    N = int(N)
    add_zero = bool(add_zero)
    symmetric_about_zero = bool(symmetric_about_zero)
    #  ====================================================

    # input check =========================================
    # vmin
    if vmin == 0.:
        raise ValueError(
                '`vmin` must be truly positive. Use `add_zero` option'
                + ' if you want to include 0.'
                )
    if not (vmin > 0.):
        raise ValueError('`vmin` must be positive.')

    # vmax
    if not vmax > vmin:
        raise ValueError('`vmax` must be greater than `vmin`.')

    # N
    if symmetric_about_zero and (not N >= 4):
        raise ValueError('`N` must be at least 4 for symmetric range.')
    elif add_zero and (not N >= 3):
        raise ValueError('`N` must be at least 3 for range with zero.')
    elif (not N >= 1):
        raise ValueError('`N` must be at least 2.')
    #  ====================================================

    if symmetric_about_zero:
        if add_zero and N % 2 == 0:
            N += 1
        elif not add_zero and N % 2 == 1:
            N += 1

        Nhalf = N // 2
        vticks_half = logarithmic_range(vmin, vmax, Nhalf+1, add_zero=True)
        if add_zero:
            vticks_pos = vticks_half
        else:
            vticks_pos = vticks_half[1:]
        vticks_neg = - vticks_half[1:][::-1]
        vticks = np.concatenate((vticks_neg, vticks_pos))
        return vticks

    if add_zero:
        # 1st point: 0
        # N-1 points: positive and logarithmic
        vticks_log = logarithmic_range(vmin, vmax, N-1, False)
        vticks = np.concatenate(([0], vticks_log))
        return vticks

    f = np.arange(N) / (N - 1)
    vticks = vmin * (vmax/vmin)**f

    return vticks
