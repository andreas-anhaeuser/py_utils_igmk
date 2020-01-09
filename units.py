def get_scale_factor(units):
    if units == 'm':
        return 1.

    if units == 'km':
        return 1000.

    raise ValueError('Unknown units: %s' % units)
