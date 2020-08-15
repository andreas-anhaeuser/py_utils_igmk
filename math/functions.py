import numpy as np

def gauss(fwhm=None, center=0., amplitude=1., sigma=None):
    """Return a Gauss function.
    
        Parameters
        ----------
        fwhm : float or None
            full width at half maximum
        sigma : float or None
            mutually exclusive with fwhm

        center : float
            (default: 0.)
        amplitude : float
            integral over the function (default: 1.)

        Returns
        -------
        callable
    """
    if sigma is None and fwhm is None:
        raise ValueError('Either fwhm or sigma must be given, but not both.')

    if sigma is not None and fwhm is not None:
        raise ValueError('Either fwhm or sigma must be given, but not both.')

    if sigma is None:
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))

    A  = amplitude / (sigma * np.sqrt(2*np.pi)) # amplitude factor
    C  = -1 / (2 * sigma**2)                    # exponent factor
    x0 = center
    f = lambda x: A * np.exp(C * (x-x0)**2)
    return f
    


def gauss_old():
    A  = amplitude / (width * np.sqrt(2*np.pi)) # amplitude factor
    C  = -1 / (2 * width**2)                    # exponent factor
    x0 = center
    f = lambda x: A * np.exp(C * (x-x0)**2)
    return f
