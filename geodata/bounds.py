"""Bounds utils."""

# PyPI modules
from shapely.geometry import Polygon
import fiona
import geopandas as gp

def convert_bounds_to_gdf(bounds, crs=None):
    """Try to convert the input to a GeoDataFrame.

        Paramters
        ---------
        bounds : anything that can be interpreted as a bounding polygon
            - lonmin, latmin, lonmax, latmax
            - a polygon
            - a GeoSeries
            - a GeoDataFrame
        crs : coordinate reference system, optional
            regular lonlat will be assumed as default

            
        Returns
        -------
        bounds : GeoDataFrame
    """
    # GeoDataFrame -> nothing to do
    # ------------------------------------------------
    if isinstance(bounds, gp.GeoDataFrame):
        return bounds

    # GeoSeries -> GeoDataFrame
    # ------------------------------------------------
    if isinstance(bounds, gp.GeoSeries):
        if crs is None:
            crs = fiona.crs.from_epsg(4326) 

        gdf = gp.GeoDataFrame(geometry=bounds)
        gdf.crs = crs

        return gdf

    # Polygon -> GeoSeries
    # ------------------------------------------------
    if isinstance(bounds, Polygon):
        geoseries = gp.GeoSeries(bounds)
        return convert_bounds_to_gdf(geoseries)

    # limits -> Polygon
    # ------------------------------------------------
    xmin, ymin, xmax, ymax = bounds
    coordinates = (
            (xmin, ymin),
            (xmax, ymin),
            (xmax, ymax),
            (xmin, ymax),
            )
    polygon = Polygon(coordinates)
    return convert_bounds_to_gdf(polygon)

def convert_bounds_to_tuple(bounds):
    if isinstance(bounds, gp.GeoDataFrame):
        return tuple(bounds.bounds.iloc[0])

    if isinstance(bounds, gp.Geoseries):
        return tuple(bounds.bounds.iloc[0])

    bounds = tuple(bounds)
    assert len(bounds) == 4
    assert bounds[2] >= bounds[0]
    assert bounds[3] >= bounds[1]

    return bounds
