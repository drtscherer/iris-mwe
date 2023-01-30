import logging
from . import logging_config
logger = logging.getLogger(__name__)
logger.setLevel(logging_config.level)

import pyproj
import math
from shapely.ops import transform
from datetime import datetime, timedelta


JDAY2000_EPOCHE = datetime(2000,1,1)

class ConversionException(Exception):
    """:meta private:"""
    def __init__(self, args):
        logger.error(args)

def convert_wgs_to_utm_zone(lon=None, lat=None, geometry=None):
    if geometry is not None:
        lon, lat = list(geometry.centroid.coords)[0]
    elif geometry is None and None in [lon,lat]:
        raise ConversionException('Either Geometry or full pair of coordinates reqired for utm conversion')
    utm_band = (math.floor((lon + 180) / 6 ) % 60) + 1
    epsg_code = '326' if lat >= 0 else '327'
    return '{}{utm:02d}'.format(epsg_code, utm=int(utm_band))

def convert_utm_geometry_to_wgs(geometry, epsg_utm):
    wgs84 = pyproj.CRS('EPSG:4326')
    utm = pyproj.CRS(f'EPSG:{epsg_utm}')
    project_utm_to_wgs84 = pyproj.Transformer.from_crs(utm, wgs84, always_xy=True).transform
    geometry_wgs = transform(project_utm_to_wgs84, geometry)
    return geometry_wgs

def convert_wgs_geometry_to_utm(geometry, epsg_utm):
    wgs84 = pyproj.CRS('EPSG:4326')
    utm = pyproj.CRS(f'EPSG:{epsg_utm}')
    project_wgs84_to_utm = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform
    geometry_utm = transform(project_wgs84_to_utm, geometry)
    return geometry_utm

def convert_jday2000_to_datetime(jday):
    return JDAY2000_EPOCHE + timedelta(days=jday)

def convert_jday_to_datetime(jday, epoche):
    return epoche + timedelta(days=jday)

def convert_jseconds_to_datetime(jseconds, epoche):
    return epoche + timedelta(seconds=jseconds)
