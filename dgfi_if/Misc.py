import logging, os
from . import logging_config
logger = logging.getLogger(__name__)
logger.setLevel(logging_config.level)

from .Conversions import convert_wgs_to_utm_zone, convert_wgs_geometry_to_utm

from typing import Union, Optional, List, Dict

import numpy as np
import pandas as pd
import geopandas as gpd
from osgeo import ogr
from datetime import datetime

import hashlib
import pyproj
import math
from shapely.ops import transform, split
from shapely.geometry import Point, LineString, MultiLineString
# import tum_design as tum

class MiscException(Exception):
    """:meta private:"""
    def __init__(self, args):
        logger.error(args)

def great_circle_dist(lon_list : Union[List[float],float], lat_list : Union[List[float],float], query_lon : float, query_lat : float) -> float:
    query_lon, query_lat = np.radians(query_lon), np.radians(query_lat)
    lon_list, lat_list = np.radians(lon_list), np.radians(lat_list)

    dlon = lon_list - query_lon
    dlat = lat_list - query_lat

    a = np.sin(dlat / 2.0)**2 + np.cos(query_lat) * np.cos(lat_list) * np.sin(dlon / 2.0)**2
    distances = 2 * np.arcsin(np.sqrt(a))
    return distances * 6371.

def get_pfaff6_basins():
    root = os.path.abspath(os.path.dirname(__file__))
    shapefile = f'{root}/pfaff_6/Hybas6.shp'
    gdf = gpd.read_file(shapefile).set_index('PFAF_ID',drop=False)
    gdf.set_crs(epsg=4326, inplace=True, allow_override=True)
    return gdf

def get_pfaff4_basins():
    root = os.path.abspath(os.path.dirname(__file__))
    shapefile = f'{root}/pfaff_4/Hybas4.shp'
    gdf = gpd.read_file(shapefile).set_index('PFAF_ID',drop=False)
    gdf.set_crs(epsg=4326, inplace=True, allow_override=True)
    return gdf

# PFAFF_6_GDF = get_pfaff6_basins()

def get_pfaff_code(lon : float, lat : float, level : int = 2, max_dist : float = 10000, raise_exception : bool = True) -> int:
    """
    Get the queried Pfaffstetter Level Code between 1 and 6 of the basin enclosing the query coords

    Parameter:
    ------
    lon : float
        longitude of POI
    lat : float
        latitude of POI
    level : int (Default 2)
        Pfaffstetter level (between 1 and 6)
    max_dist : float (Default 10000)
        maximum distance [m] from basin. Used when given coords are not within a pfaffstetter level 6 basin. E.g estuaries.
    raise_exception : bool
        If True raises an exception when no basin was found, else None is returned

    Returns
    ------
    int
        Pfaffstetter Level X Code
    """
    if level not in [1,2,3,4,5,6]:
        raise MiscException(f'Pfaffstetter Level {level} not supported.')
    poi = Point(lon,lat)
    gdf = PFAFF_6_GDF[PFAFF_6_GDF.contains(poi)]
    pfaf_id =None
    if gdf.empty:
        distances = PFAFF_6_GDF.distance(poi)
        max_degree_dist = 0.00001 * max_dist
        distances = distances[distances <= max_degree_dist]
        if len(distances.index) > 0:
            epsg_code = convert_wgs_to_utm_zone(lon,lat)
            utm_poi = convert_wgs_geometry_to_utm(poi, epsg_code)
            gdf = PFAFF_6_GDF[PFAFF_6_GDF.index.isin(distances.index)]
            gdf.to_crs(epsg=epsg_code,inplace=True)
            distances = gdf.distance(utm_poi).sort_values()
            if distances.iloc[0] <= max_dist:
                pfaf_id = PFAFF_6_GDF.PFAF_ID.loc[distances.index[0]]
                logger.warning(f'Coordinates {lon} {lat} not within Pfaffstetter Level 6 Basins. Returned nearest Basin. Distance: {distances.iloc[0] / 1000:.2f} km')
    else:
        pfaf_id = gdf.PFAF_ID.iloc[0]
    if pfaf_id is None and raise_exception:
        raise MiscException(f'Coordinates {lon} {lat} not within Pfaffstetter Level 6 Basins.')
    elif pfaf_id is None:
        return None
    else:
        return int(pfaf_id / (10 ** (6-level)))

def get_pfaff2_code(lon : float, lat : float) -> int:
    """
    Get the Pfaffstetter Level 2 Code of the basin enclosing the query coords

    Parameter:
    ------
    lon : float
        longitude of POI
    lat : float
        latitude of POI

    Returns
    ------
    int
        Pfaffstetter Level 2 Code
    """
    root = os.path.abspath(os.path.dirname(__file__))
    shapefile = f'{root}/pfaff_2/pfaff_2.shp'
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapefile, 0)
    layer = dataSource.GetLayer()
    pfaff_2_code = None

    wkt = f'POINT ({lon} {lat})'
    layer.SetSpatialFilter(ogr.CreateGeometryFromWkt(wkt))

    for feature in layer:
        pfaff_2_code = feature.GetField("PFAF_ID")
        break
    if pfaff_2_code == None:
        pfaff_2_code = get_pfaff_code(lon,lat,2)

    driver, dataSource, layer = None, None, None
    return pfaff_2_code


def md5hash_str(string):
    return hashlib.md5(string.encode()).hexdigest()

def spherical_distance(lon1 : Union[float, np.ndarray], lat1 : Union[float, np.ndarray], lon2 : Union[float, np.ndarray], lat2 : Union[float, np.ndarray], ellipsoid='wgs84') -> Union[float, np.ndarray]:
    ellipsoids = {}
    ellipsoids['wgs84'] = 6378.1370	
    ellipsoids['grs80'] = 6378.1370
    ellipsoids['topex'] = 6378.1363

    R = ellipsoids[ellipsoid]
    dLat = np.radians(lat2-lat1)
    dLon = np.radians(lon2-lon1)
    a = np.sin(dLat/2) * np.sin(dLat/2) +np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dLon/2) * np.sin(dLon/2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return R * c
    
def timeseries_to_dataframe(timeseries):
    dates = []
    for _ts in timeseries.values():
        if _ts is not None:
            dates += list(_ts.keys())
    dates = sorted(set(dates))
    for _ts_key, _ts in timeseries.items():
        if _ts is not None:
            timeseries[_ts_key] = [_ts[date] if date in _ts.keys() else np.nan for date in dates]
        else:
            timeseries[_ts_key] = [np.nan] * len(dates)
    dates_str = dates
    dates = [datetime.strptime(d, "%Y-%m-%d") for d in dates_str]
    data = {'dates': dates, 'dates_str':dates_str}
    data.update(timeseries)
    return pd.DataFrame(data).set_index('dates')

def validate_geometry(geometry):
    if geometry.geom_type not in ['LineString', 'MultiLineString','Polygon']:
        raise MiscException('Invalid Geometry!')
    if not geometry.is_valid:
        if geometry.geom_type == 'LineString':
            lon, lat = zip(*geometry.coords)
            lon, lat = lon[0::20], lat[0::20]
            geometry = LineString(list(zip(lon,lat)))
        elif geometry.geom_type == 'MultiLineString':
            coords = []
            for geom in geometry.geoms:
                lon, lat = zip(*geom.coords)
                lon, lat = lon[0::20], lat[0::20]
                coords.append(list(zip(lon,lat)))
            geometry = MultiLineString(coords)
        elif geometry.geom_type == 'Polygon':
            geometry = geometry.buffer(0)
    if not geometry.is_valid:
        raise MiscException('Invalid Geometry!')
    return geometry

def validate_line(geometry):
    logger.warning('validate_line is depricated, use validate_geometry!')
    return validate_geometry(geometry)

def get_geom_coordinates(geometry):
    "Returns list of Coordinates in Shapely Objects"
    if geometry.geom_type not in ['LineString', 'MultiLineString','MultiPoint','Polygon','MultiPolygon']:
        raise MiscException(f'Invalid Geometry! ({geometry.geom_type})')
    if geometry.geom_type == 'LineString':
        coords = list(geometry.coords)
    elif geometry.geom_type == 'Polygon':
        coords = list(geometry.exterior.coords)
    elif geometry.geom_type == 'MultiPolygon':
        coords = []
        for geom in geometry.geoms:
            coords += list(geom.exterior.coords)
    elif geometry.geom_type in  ['MultiLineString','MultiPoint']:
        coords = []
        for geom in geometry.geoms:
            coords += list(geom.coords)
    return coords

def linestring_helper(lon : List[float], lat : List[float], wkt : bool = True) -> Union[str,LineString]:
    """
    Convert lists of lon and lat to Linestring WKT.
    Converts lon to -180 to 180 if necessary.
    Splits Lines at Anti-Meridian.
    """
    if len(lon) < 2:
        raise MiscException('There must be more than one Coord to create a Line')
    lon = [x - 360. if x > 180 else x for x in lon]
    if lon[0] > lon[-1]:
        lon = [(x + 360.) % 360 for x in lon]
        geometry = LineString(list(zip(lon,lat)))
        splitted = split(geometry,LineString([[180,-90],[180,90]]))
        coords = []
        for geom in splitted.geoms:
            _lon, _lat = zip(*geom.coords)
            if max(_lon) > 180:
                _lon = [x - 360. for x in _lon]
            coords.append(list(zip(_lon,_lat)))
        geometry = MultiLineString(coords)
    else:
        geometry = LineString(list(zip(lon,lat)))
    if wkt:
        return geometry.wkt
    else:
        return geometry

def icesat2_along_get_conf_threshold_by_angle(angle,max_confidence,max_angle):
    if angle > 90:
        angle = 180-angle
    if angle < max_angle:
        return max_confidence-(max_confidence/max_angle)*angle
    else:
        return 0

def openadb_spherical_distance(lon1, lat1, lon2, lat2, ellipsoid='wgs84'):
		
    ' mean equator radius '
    ellipsoids = {}
    ellipsoids['wgs84'] = 6378.1370	
    ellipsoids['grs80'] = 6378.1370
    ellipsoids['topex'] = 6378.1363

    R = ellipsoids[ellipsoid];
    dLat = math.radians(lat2-lat1);
    dLon = math.radians(lon2-lon1);
    a = math.sin(dLat/2) * math.sin(dLat/2) +math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dLon/2) * math.sin(dLon/2); 
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)); 
    d = R * c;

    # distance in km!"
    return d