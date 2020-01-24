# standard modules
import collections

# PyPI modules
import numpy as np

def round_digits(x, digits=1):
    if digits <= 0:
        raise ValueError('`digits` must be a positive integer.')

    if isinstance(x, collections.Iterable):
        N = len(x)
        return np.array([round_digits(x[n], digits) for n in range(N)])

    if x == 0:
        return x

    magnitude_fractional = np.log10(np.abs(x))
    magnitude = int(np.floor(magnitude_fractional))
    decimals = digits - 1 - magnitude
    return np.round(x, decimals)
