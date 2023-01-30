import logging
from . import logging_config
logger = logging.getLogger(__name__)
logger.setLevel(logging_config.level)
from . import Utilities
from OpenADB.MVA import read_MVA, mask_MVA, get_index_by_polygon
from OpenADB.Time import julianDayDate

import os, re, sys
from typing import Optional, List, Union
import h5py
import numpy as np
import pandas as pd
import sqlite3
from multiprocessing import Pool
import multiprocessing as mp
import math
import pathlib
from shapely.geometry import MultiLineString, LineString, Polygon
import shapely
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map
from datetime import datetime

DEFAULT_MVA_ROOT = ...
DEFAULT_CACHE_ROOT = ...
DEFAULT_ORIGINAL_ICESAT2_DATA_ROOT = ...
ICESAT2_EPOCHE = datetime(2018,1,1)

class MVAException(Exception):
    """:meta private:"""
    def __init__(self, args):
        logger.debug(args)

class MVA_IF:
    def __init__(self,use_cache=False,cache_root=None,mva_root=None):
        r"""
        Convenience interface class for interactions with the DGFI MVA data.

        Parameters
        ----------
        use_cache : boolean
            set True if you want to use cached data if available
        cache_root : string
            cache_root of the cached orbit data
        mva_root : string
            path to the mva root
        """
        self.cache_root = pathlib.Path(cache_root) if cache_root is not None else DEFAULT_CACHE_ROOT
        self.mva_root = pathlib.Path(mva_root) if mva_root is not None else DEFAULT_MVA_ROOT
        self.use_cache = use_cache
        if not self.mva_root.exists():
            raise MVAException(f'Given MVA root ({str(self.mva_root)}) does not exist')

    def init_cache(self, mission):
        r"""
        Debug only
        """
        ...

    def delete_cache(self, mission=None):
        r"""
        Debug only
        """
        ...

    def check_cache(self, hash, mission):
        r"""
        Debug only
        """
        ...

    def read_cache(self, hash, mission):
        r"""
        Debug only
        """
        ...

    def insert_cache(self,hash,mission,pass_nr,cycle_nr,geometry):
        r"""
        Debug only
        """
        ...

    def insert_cache_many(self, input_list, mission):
        r"""
        Debug only
        """
        ...

    def get_orbit(
        self,
        mission : str,
        pass_nr : str,
        cycle_nr : str,
        read_cache : Optional[bool] = None,
        write_cache : Optional[bool] = None,
        lat_bounds : Optional[tuple] = None,
        ):
        r"""
        Convenience function to extract the orbit geometry of a mission, path, cylce from MVA in standard longitude format (-180 to 180)
        
        Parameters
        ----------
        mission : string
            The mva mission identifier eg. 'jason2_hf'
        pass_nr : string
            The mva mission pass number with correct leading zeroes eg. '0042'
        cycle_nr : string
            The mva mission cycle number with correct leading zeroes eg. '042'
        read_cache : boolean
            If True, checks the cache for already extracted data. Default value is taken from the class use_cache parameter.
        write_cache : boolean
            If True, writes the extracted data to cache. Default value is taken from the class use_cache parameter.

        RETURNS
        ----------
        WKT string
            Well Known Text representation of the orbit geometry. Linestring if anti meridian is not crossed. Otherwise splitted into MultiLineString at anti meridian. Longitude format is -180 to 180.
            Returned geometry may be invalid but contains all MVA hf points.
        """
        
        read_cache = read_cache if read_cache is not None else self.use_cache
        write_cache = write_cache if write_cache is not None else self.use_cache
        in_cache = False
        if any([read_cache,write_cache]):
            cache_location = self.cache_root.joinpath(f'{mission}_cache.db')
            if not cache_location.exists():
                self.init_cache(mission)
            hash = Utilities.md5hash_str(f'{mission}{pass_nr}{cycle_nr}')
            in_cache = self.check_cache(hash,mission)
            
        
        if read_cache is False or in_cache is False:
            records = None
            fn = f'{cycle_nr}_{pass_nr}orbit.00'
            orbit_path = str(self.mva_root / mission / cycle_nr / fn)
            if lat_bounds is not None:
                options = {}
                options['glat'] = ('+',lat_bounds[0],lat_bounds[1])
                records = mask_MVA(orbit_path,options)

            options = {}
            options['parameters'] = ['glon','glat']
            if records is not None:
                options['records'] = records
            orbit_object = read_MVA(str(orbit_path),options)
            
            if orbit_object['num_records'] == 0:
                raise MVAException('Orbit does not exist!')

            lon, lat = [x for x in orbit_object['glon'] if not math.isnan(x)], [x for x in orbit_object['glat'] if not math.isnan(x)]
            
            if len(lon) < 2:
                raise MVAException('Orbit does not exist!')

            geometry = Utilities.linestring_helper(lon,lat)
            if write_cache:
                self.insert_cache(hash,mission,pass_nr,cycle_nr,geometry)
        else:
            geometry = self.read_cache(hash,mission)
        return geometry

    def get_orbit_parallel_helper(self, fn):
        r"""
        Internal helper function to process the results of a multiprocessing pool.
        """
        cycle_nr, pass_nr = str(fn).split('/')[-1][:-8].split('_')
        mission = str(fn).split('/')[4]
        try:
            geometry = self.get_orbit(mission,pass_nr,cycle_nr,write_cache=False)
        except MVAException:
            logger.debug(f'Error during MVA Extraction of {mission} Cycle {cycle_nr} Pass {pass_nr}. Skipping.')
            return None
        return {'cycle_nr':cycle_nr,'pass_nr':pass_nr,'mission':mission,'geometry':geometry}

    def get_passes(self,mission,pass_filter,cycle_filter=None,num_cpus=None,remaining_cpus = 1,read_cache=None,write_cache=None,silent=True):
        r"""
        Convenience function to extract multiple orbit pass geometries with metadata of a mission from MVA in standard longitude format (-180 to 180)
        
        Parameters
        ----------
        mission : string
            The mva mission identifier eg. 'jason2_hf'
        pass_filter : List of strings
            The mva mission pass numbers with correct leading zeroes eg. ['0042','0043',...]
        cycle_filter : List of strings
            The mva mission cycle numbers with correct leading zeroes eg. ['042','043',...]. Default: None (all cycles will be extracted).
        num_cpus : Integer
            Number of CPUs used in multiprocessing pool. Default: System CPUs - 1
        remaining_cpus : Integer
            Number of remaining CPUs not to be used by the multiprocessing pool. Default: 1
        read_cache : boolean
            If True, checks the cache for already extracted data. Default value is taken from the class use_cache parameter.
        write_cache : boolean
            If True, writes the extracted data to cache. Default value is taken from the class use_cache parameter.

        RETURNS
        ----------
        List of Dictionaries per pass
            cycle_nr : Cylce number of the orbit geometry
            pass_nr : Pass number of the orbit geometry
            mission : Mission number of the orbit geometry
            geometry : Well Known Text representation of the orbit geometry. Linestring if anti meridian is not crossed. Otherwise splitted into MultiLineString at anti meridian. Longitude format is -180 to 180.
            returned features geometry may be invalid but contains all MVA hf points.
        """
        read_cache = read_cache if read_cache is not None else self.use_cache
        write_cache = write_cache if write_cache is not None else self.use_cache
        if any([read_cache,write_cache]):
            cache_location = self.cache_root.joinpath(f'{mission}_cache.db')
            if not cache_location.exists():
                self.init_cache(mission)
        
        logger.debug(f'Gathering selected {mission} passes from MVA.')
        features = []
        
        cycles = [c for c in os.listdir(self.mva_root.joinpath(mission)) if re.match(r"^[0-9]+$",c)]
        files = []
        for cycle_nr in tqdm(cycles, desc=f'Searching {mission} Files', leave=False):
            if cycle_filter is not None and cycle_nr not in cycle_filter:
                continue
            for pass_nr in pass_filter:
                fn = f'{cycle_nr}_{pass_nr}orbit.00'
                p = self.mva_root.joinpath(mission).joinpath(cycle_nr).joinpath(fn)
                if p.exists():
                    files.append(p)

        features = process_map(self.get_orbit_parallel_helper, files, desc=f'Extracting {mission} Files', leave=False)
        
        features = [f for f in features if f is not None]
        insert_list = []
        if write_cache:
            for f in tqdm(features, desc=f'Preparing {mission} Files', leave=False):
                mission = f['mission']
                cycle_nr = f['cycle_nr']
                pass_nr = f['pass_nr']
                geometry = f['geometry']
                hash = Utilities.md5hash_str(f'{mission}{pass_nr}{cycle_nr}')
                insert_list.append([hash,mission,pass_nr,cycle_nr,geometry])
            if len(insert_list) > 0:
                if not silent:
                    print('', end="\r")
                    print('Writing MVA Data to Cache', end="\r")
                self.insert_cache_many(insert_list,mission)
        return features

    def get_passes_by_polygon(self, mission : str, polygon : Polygon):
        """
        Returns Passes based on mission and AOI polygon.

        Parameters
        ----------
        mission : string
            The mva mission identifier eg. 'jason2_hf'
        polygon : shapely.geometry.Polygon
            The AOI geometry polygon used to query the MVA data

        Returns:
        ----------
        List of Dictionaries per pass
            cycle_nr : Cylce number of the orbit geometry
            pass_nr : Pass number of the orbit geometry
            mission : Mission number of the orbit geometry
            geometry : Well Known Text representation of the orbit geometry. Linestring if anti meridian is not crossed. Otherwise splitted into MultiLineString at anti meridian. Longitude format is -180 to 180.
            returned features geometry may be invalid but contains all MVA hf points.
        """
        if not isinstance(polygon, Polygon):
            raise MVAException("Error in get_passes_by_polygon: Polygon must be instance of shapely.geometry.Polygon")
        if not polygon.is_valid:
            raise MVAException("Polygon must be valid")
        query_polygon = []
        for x, y in polygon.exterior.coords:
            if x < 0:
                x+=360.0
            query_polygon.append([x,y])

        # logger.info(f'Requesting {mission} data in polygon')
        result = get_index_by_polygon(mission,query_polygon,version='auto')
        
        # logger.info('Iterating Features')
        features = []
        for cycle_nr, passes in result.items():
            for pass_nr, records in passes.items():
                if len(records) < 2:
                    continue
                records = pd.DataFrame.from_dict(records,orient='index')
                geometry = Utilities.linestring_helper(records.glon, records.glat)
                features.append({'cycle_nr':cycle_nr, 'pass_nr':pass_nr, 'mission':mission, 'geometry':geometry, "records": records.index.tolist()})
        # logger.info(f'Done, {len(features)} features')
        return features

    def get_icesat2_data_by_polygon(self, polygon : Polygon, missions : Optional[Union[List[str],str]] = None, original_data_root : Optional[str] = None, silent=False) -> pd.DataFrame:
        if not isinstance(polygon, Polygon):
            raise MVAException("Error in get_icesat2_data_by_polygon: Polygon must be instance of shapely.geometry.Polygon")
        if not polygon.is_valid:
            raise MVAException("Polygon must be valid")
        result = {}
        if missions is None:
            missions = [f'icesat2_gt{x}_atl13v5_hf' for x in ['1l','1r','2l','2r','3l','3r']]
        elif isinstance(missions,str):
            missions = [missions]
        for mission in tqdm(missions, desc='Extracting ICESat-2 Missions from MVA', leave = False):
            result[mission] = []
            _, beam, dataset, _ = mission.split('_')
            features = self.get_passes_by_polygon(mission, polygon)
            for feature in tqdm(features, desc=f'MVA: Iterating {mission} pass features in Polygon', leave=False, disable=silent):
                cycle_nr, pass_nr = feature['cycle_nr'], feature['pass_nr']
                options = {}
                options['records'] = feature['records']
                fn = f'{cycle_nr}_{pass_nr}water_body.00'
                water_object = read_MVA(str(self.mva_root / mission / cycle_nr / fn),options)
                fn = f'{cycle_nr}_{pass_nr}elev.00'
                elev_object = read_MVA(str(self.mva_root / mission / cycle_nr / fn),options)
                fn = f'{cycle_nr}_{pass_nr}time.00'
                time_object = read_MVA(str(self.mva_root / mission / cycle_nr / fn),options)
                fn = f'{cycle_nr}_{pass_nr}geoh.07'
                geoid07_object = read_MVA(str(self.mva_root / mission / cycle_nr / fn),options)
                # fn = f'{cycle_nr}_{pass_nr}geoh.10'
                # geoid10_object = read_MVA(str(self.mva_root / mission / cycle_nr / fn),options)
                fn = f'{cycle_nr}_{pass_nr}cloud_flag.00'
                cloud_flag_object = read_MVA(str(self.mva_root / mission / cycle_nr / fn),options)
                fn = f'{cycle_nr}_{pass_nr}snow_ice_atl09.00'
                snow_ice_atl09_object = read_MVA(str(self.mva_root / mission / cycle_nr / fn),options)
                options['parameters'] = ['glon','glat']
                fn = f'{cycle_nr}_{pass_nr}orbit.00'
                orbit_object = read_MVA(str(self.mva_root / mission / cycle_nr / fn),options)

                jday = np.array(time_object['jday'])

                # logger.info('MVA Done')
                # strong = self.icesat2_param_helper(mission, cycle_nr, pass_nr, jday, original_data_root,silent=silent)
                strong = None
                # logger.info('Orient Done')
                
                lon, lat = np.array(orbit_object['glon']), np.array(orbit_object['glat'])
                lon = np.where(lon > 180, lon - 360, lon)
                alongtrack_distance = Utilities.spherical_distance(lon[0],lat[0],lon,lat)

                feature['mission'] = mission
                feature['beam'] = beam
                feature['data'] = pd.DataFrame({
                    'lon': lon,
                    'lat': lat,
                    'wgs_coord': list(zip(lon,lat)),
                    'alongtrack_distance' : alongtrack_distance * 1000.,
                    'elev': np.array(elev_object['elev']),
                    'depth': np.array(water_object['water_depth']),
                    'atl13_slope': np.array(water_object['segment_slope_trk_bdy']),
                    'jday': jday,
                    'strong': strong,
                    'geoh07': np.array(geoid07_object['geoh']),
                    # 'geoh10': np.array(geoid10_object['geoh']),
                    'water_body_id': np.array(water_object['inland_water_body_id']),
                    'water_body_type': np.array(water_object['inland_water_body_type']),
                    'cloud_flag_asr_atl09': np.array(cloud_flag_object['cloud_flag_asr_atl09']),
                    'cloud_flag_atm_atl09': np.array(cloud_flag_object['cloud_flag_atm_atl09']),
                    'layer_flag_atl09': np.array(cloud_flag_object['layer_flag_atl09']),
                    'qf_cloud': np.array(cloud_flag_object['qf_cloud']),
                    'qf_ice': np.array(cloud_flag_object['qf_ice']),
                    'snow_ice_atl09': np.array(snow_ice_atl09_object['snow_ice_atl09']),
                    })
                result[mission].append(feature)
        return result