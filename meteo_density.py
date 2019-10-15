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
