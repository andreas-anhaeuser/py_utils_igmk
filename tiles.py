#!/usr/bin/python
"""Tiles of 2d rasters."""

# standard modules
import itertools

# PyPI modules
import numpy as np

class Tile(object):
    """A 2d tile."""
    def __init__(self, ymin, ymax, xmin, xmax):
        self.ymin = ymin
        self.ymax = ymax
        self.xmin = xmin
        self.xmax = xmax

    def __str__(self):
        """Return str representation."""
        return (
                'Tile with y: [%i, %i), x: [%i, %i)'
                % (self.ymin, self.ymax, self.xmin, self.xmax)
                )

    def __repr__(self):
        """Return a str."""
        return str(self)

    def shape(self):
        """Return a tuple of int."""
        Ny = self.ymax - self.ymin
        Nx = self.xmax - self.xmin
        return (Ny, Nx)

    def size(self):
        """Return an int."""
        shape = self.shape()
        return np.prod(shape)

    def zeros(self, dtype=float):
        shape = self.shape()
        return np.zeros(shape, dtype=dtype)

    def ones(self, dtype=float):
        shape = self.shape()
        return np.ones(shape, dtype=dtype)

    def nans(self):
        shape = self.shape()
        return np.nan * np.ones(shape)

    def get_index(self):
        idx_y = slice(self.ymin, self.ymax)
        idx_x = slice(self.xmin, self.xmax)
        return (idx_y, idx_x)

def get_tile_numbers(edge_size_y, edge_size_x, total_tile_number):
    total_area = edge_size_y * edge_size_x
    area_per_tile = total_area / total_tile_number
    edge_size_per_tile = np.sqrt(area_per_tile)
    Ntiles_y = np.round(edge_size_y / edge_size_per_tile)
    Ntiles_x = np.round(total_tile_number / Ntiles_y)
    return int(Ntiles_y), int(Ntiles_x)

def tile_2d_raster(Ny, Nx, max_pixels_per_tile):
    """Return list of tiles that are not large than `pixels_per_tiles.
    
        Parameters
        ----------
        Ny : int
            raster size in y-direction
        Nx : int
            raster size in x-direction
        max_pixels_per_tile : int
            maximum number of pixels per tile

        Returns
        -------
        tiles : list of Tile
    """
    if Ny * Nx <= max_pixels_per_tile:
        # return list of one Tile spanning the whole domain
        return [Tile(0, Ny, 0, Nx)]

    ceil = lambda x : int(np.ceil(x))
    ppt = int(max_pixels_per_tile)

    # y-size of tile
    edge_size = np.sqrt(ppt)
    Ntiles_y = ceil(Ny / edge_size)
    tile_size_y = ceil(Ny / Ntiles_y)

    # x-size of tile
    N_per_y_line = tile_size_y * Nx
    Ntiles_x = ceil(N_per_y_line / ppt)
    tile_size_x = ceil(Nx / Ntiles_x)
    
    tiles = []
    it = itertools.product(range(Ntiles_y), range(Ntiles_x))
    for ntile_y, ntile_x in it:
        ymin = ntile_y * tile_size_y
        ymax_naive = ymin + tile_size_y
        ymax = min(ymax_naive, Ny)

        xmin = ntile_x * tile_size_x
        xmax_naive = xmin + tile_size_x
        xmax = min(xmax_naive, Nx)

        tile = Tile(ymin, ymax, xmin, xmax)
        tiles.append(tile)

    return tiles
