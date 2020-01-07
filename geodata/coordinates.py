import numpy as np

def compute_coordinates_from_transform(meta):
    """Return x and y as two 1d-arrays."""
    Nx = meta['width']
    Ny = meta['height']
    transform = meta['transform']

    # i : indices in x-direction
    # j : inideces in y-direction
    i = np.arange(Nx)
    j = np.arange(Ny)

    # apply affine transformation from index to coordinate
    # `x0` and `y0` remain unused
    # (the zeros can be replace by any value without altering the output)
    x, y0 = transform * (i, 0)
    x0, y = transform * (0, j)
    return x, y
