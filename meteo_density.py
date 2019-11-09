import numpy as np

_dobson_unit = 2.2687e20            # m-2
_avogadro_constant = 6.02214076e23  # mol-1

def volume_to_mass_density(v, M, p=1013.25e2, T=288.15):
    """

        Parameters
        ----------
        v : float
            (m3 m-3) volume density
        M : float
            (kg) molecular mass
        p : float
            (N m-2) pressure
        T : float
            (K) temperature

        Returns
        -------
        rho : float
            (kg m-3) mass density
    """
    kB = 1.380649e-23   # (J K-1) Boltzmann constant
    rho = p * M / (T * kB) * v
    return rho

def si_to_dobson_unit(value, substance):
    """Convert from kg m-2 to Dobson unit."""
    m = molecular_mass(substance)
    N_molec = value / m             # (m-2) molecules per surface area
    return N_molec / _dobson_unit

def molecular_mass(substance):
    if substance.lower() == 'no2':
        M = 46.006e-3
    elif substance.lower() == 'so2':
        M = 64.066e-3
    else:
        raise NotImplementedError()
    return M / _avogadro_constant
