# PyPI modules
import numpy as np

# constants
_avogadro_constant = 6.02214076e23  # mol-1
_boltzmann_constant = 1.380649e-23  # (J K-1)
_dobson_unit = 2.2687e20            # m-2

# defaults
_pressure = 1013.25e2   # N m-2
_temperature = 288.15   # K

def ppb_to_ug(ppb, substance, p=None, T=None):
    M = molecular_mass(substance)
    v = ppb * 1e-9
    mass_si = volume_to_mass_density(v, M, p, T)
    mass_ug = mass_si * 1e9
    return mass_ug

def volume_to_mass_density(v, M, p=None, T=None):
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
    # defaults ============================================
    if p is None:
        p = _pressure
    if T is None:
        T = _temperature
    # =====================================================

    rho = p * M / (T * _boltzmann_constant) * v

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
