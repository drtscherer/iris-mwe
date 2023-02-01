import logging
from . import logging_config
logger = logging.getLogger(__name__)
logger.setLevel(logging_config.level)

import pyproj
import math
from shapely.ops import transform
from datetime import datetime, timedelta
import calendar

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

def julianDayDate(jd2000):
    """Function for the conversion from Julian Day 2000 to date

    @param jd2000 (float) - Julian Day 2000
    @return object - Date object
        
    """    

    cumday = [0,31,59,90,120,151,181,212,243,273,304,334,365];
    t = jd2000 + 0.5
    ' tmp is created because of float uncertainties '
    #~ print (t,int(t))
    tmp = float('%.8f' % ((t -int(t))))
    sec = tmp*86400.

    if int(sec) == 86400: # fix because of rounding when jday == 1000.499999999993
        sec = 0
        t += 1

    #~ print (sec)
    #	print t,sec 
    if sec < 0.:
        sec = sec + 86400.
    if t >= 0.:
        n4 = int(t/1461)
    else:
        n4 = int(t/1461) - 1
    #	print n4,sec,t
    if t >= 0.:
        it = 1
    else:
        it = 0
    it = it + int(t) - n4*1461;
    #	print it
    if it == 60:
        day = 29
        m = 2
        year = 2000 + 4*n4
    else:
        if it > 60:
            it = it -1
        n1 = int((it -1) /365)
        it = it -n1*365
    #	print n1,it
        m = int((it -1)/31)
        while it > cumday[m]:
            m = m + 1
        day = it - cumday[m - 1]
    #	print n4,n1
        year = 2000 + 4*n4 +n1
    hour = int(sec/3600)
    sec = sec - hour*3600
    minute = int(sec/60)	
    sec = sec - minute*60	
    month = m

        
    if calendar.isleap(year):		
        days = 0
        cumday = [0,31,60,91,121,152,182,213,244,274,305,335,366];
        days =  cumday[month-1]
        days = days + day - 1
        dezDate = (days*24*60*60+hour*60*60+minute*60+sec)/(366*24*60*60)+year
    else:		
        days = 0
        cumday = [0,31,59,90,120,151,181,212,243,273,304,334,365];
        days =  cumday[month-1]
        days = days + day - 1		
        dezDate = (days*24*60*60+hour*60*60+minute*60+sec)/(365*24*60*60)+year		

    object = {}
    object['year'] = year
    object['month'] = month
    object['day'] = day
    object['hour'] = hour
    object['minute'] = minute
    object['second'] = sec
    object['doy'] = datetime(year,month,day,0,0,0).strftime('%j')
    object['date'] = "%04d-%02d-%02d %02d:%02d:%02d" % (year,month,day,hour,minute,sec)
    object['date_short'] = "%04d%02d%02d%02d%02d%02d" % (year,month,day,hour,minute,sec)
    return object