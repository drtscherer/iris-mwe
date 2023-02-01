import logging
from . import logging_config
logger = logging.getLogger(__name__)
logger.setLevel(logging_config.level)
import pathlib

import h5py
import numpy as np
import pandas as pd
import geopandas as gpd
# import contextily as cx
from osgeo import ogr
import fiona, fiona.crs
from . import Utilities, Conversions
# from .MVA_IF import MVAException
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR, LinearSVR
# from OpenADB.Tools import spherical_distance
# from OpenADB.Time import julianDayDate
import datetime
import pickle

from shapely.geometry import MultiLineString, LineString, Point, Polygon, mapping, box
from shapely import wkt
from shapely.ops import unary_union, nearest_points
from typing import Dict, List, Optional, Union
from tqdm.auto import tqdm
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
# from matplotlib_scalebar.scalebar import ScaleBar
import warnings
from copy import deepcopy

# import tum_design as tum

class SWORDException(Exception):
    """:meta private:"""
    def __init__(self, args):
        logger.debug(args)
class SWORDConnectionException(Exception):
    """:meta private:"""
    def __init__(self, args):
        logger.debug(args)

class SWORD_IF:  
    def __init__(self, sword_root=None,sword_version=None):
        self.sword_version = sword_version if sword_version is not None else '2.0'
        self.sword_base_root = pathlib.Path(sword_root) if sword_root is not None else ...
        self.sword_root = self.sword_base_root.joinpath(f'v{self.sword_version}/Reaches_Nodes/')
        self.dset = None
        self.pfaff_2_code = None
        self.continent = None
        self.nc_path = None
        self.reach_id = None
        if not self.sword_root.exists():
            raise SWORDException(f'Given SWORD root ({str(self.sword_root)}) does not exist')

    def list_reaches(self):
        dic = {}
        root = self.sword_root.joinpath('netcdf')
        for file in root.iterdir():
            if '.nc' in file.name:
                dset = self.get_dset(file)
                dic[file.name] = dset['reaches']['reach_id'][:]
        return dic

    def get_pfaff2_code(self, lon, lat) -> int:
        """Returns Pfaffstetter Level 2 basin at given coord required to locate correct SWORD file."""
        pfaff_2_code = Utilities.get_pfaff_code(lon,lat,2)
        return pfaff_2_code

    def search_continent(self, pfaff_2_code) -> str:
        """ Returns Continent of given Pfaffstetter Level 2 basin required to locate correct SWORD file"""
        root = self.sword_root.joinpath('shp')
        continent = None
        for path in root.iterdir():
            if not path.is_dir():
                continue
            for file in path.iterdir():
                if f'hb{pfaff_2_code}' in file.name:
                    continent = path.name
                    break
        if continent is None:
            raise SWORDException(f'Continent not found. Invalid Pfaff 2 Doce? ({pfaff_2_code})')
        return continent

    def get_nodes_in_aoi(self, aoi: Polygon, parameters=['node_id','reach_id']) -> Dict[str, List]:
        r"""
        Returns all required parameters of all SWORD nodes within given aoi

        Parameters
        ----------
        aoi : shapely.geometry.Polygon
            Area of Interest
        parameters : list of strings (Default: ['node_id','reach_id'])
            returned paramters of the nodes within the aoi.

        Returns
        ----------
        dictionary
            queried parameters : values within aoi
        """
        poi = aoi.centroid
        minx, miny, maxx, maxy = aoi.bounds
        nc_path = self.get_dset_path(lon=poi.x,lat=poi.y)
        dset = self.get_dset(nc_path)

        inside = np.all([dset['nodes']['x'][:] >= minx, dset['nodes']['y'][:] >= miny, dset['nodes']['x'][:] <= maxx, dset['nodes']['y'][:] <= maxy],axis=0)
        for i in np.where(inside)[0]:
            if not aoi.contains(Point(dset['nodes']['x'][i],dset['nodes']['y'][i])):
                inside[i] = False
        nodes = {}
        for parameter in parameters:
            if parameter not in dset['nodes'].keys():
                raise SWORDException(f'Queried paramter "{parameter}" not in dset. Available parameters: {(", ").join(dset["nodes"].keys())}')
            nodes[parameter] = dset['nodes'][parameter][:][inside].tolist()
        dset.close()
        dset = None
        return nodes

    def get_nodes_in_reach(self, reach_id : int, parameters=['node_id','x','y']) -> Dict[str, List]:
        r"""
        Returns all required parameters of all SWORD nodes within given reach

        Parameters
        ----------
        reach_id : int
            reach id of query reach
        parameters : list of strings (Default: ['node_id','x','y')
            returned paramters of the nodes within the aoi.

        Returns
        ----------
        dictionary
            queried parameters : values within aoi
        """
        nodes = {}
        pfaff_2_code = int(str(reach_id)[:2])
        nc_path = self.get_dset_path(pfaff_2_code=pfaff_2_code)
        self.dset = self.get_dset(nc_path)
        mask = self.dset['nodes']['reach_id'][:] == int(reach_id)
        for parameter in parameters:
            if parameter not in self.dset['nodes'].keys():
                raise SWORDException(f'Queried paramter "{parameter}" not in dset. Available parameters: {(", ").join(dself.set["nodes"].keys())}')
            nodes[parameter] = self.dset['nodes'][parameter][:][mask].tolist()
        self.close_dset()
        return nodes


    def search_reach_id(self, lon : float, lat : float, decision_criteria : str = 'facc', return_node = False, catch_exceptions = False, max_dist : float = np.inf) -> int:
        r"""
        Returns closest reach id to given coordinates

        Parameters
        ----------
        lon : float
            Longitude of query location
        lat : float
            latitude of query location
        decision_criteria : str (Default 'facc'')
            Parameter of SWORD reach data set to decide which reach should be returned in case multiple reaches have the same distance.
            By default the reach wchich the largest estimated flow will be returned.
        max_dist : float (Default np.inf)
            maximum distance to centerline in km, if further returns None
        Returns
        ----------
        reach_id : int
        """
        reach_id = None
        node_id = None
        decision_value = 0
        try:
            nc_path = self.get_dset_path(lon=lon,lat=lat)
            dset = self.get_dset(nc_path)

            # calc great circle distance
            nodes_lon, nodes_lat = dset['nodes']['x'][:], dset['nodes']['y'][:]
            distances = Utilities.great_circle_dist(nodes_lon, nodes_lat, lon, lat)
            if np.amin(distances) > max_dist:
                return None
            closest = np.where(distances == np.amin(distances))
            if len(closest) == 1:
                reach_id = dset['nodes']['reach_id'][closest][0]
                node_id = dset['nodes']['node_id'][closest][0]
            else:
                for i in closest:
                    if _val := dset['nodes'][decision_criteria][i] > decision_value:
                        decision_value = _val
                        reach_id = dset['nodes']['reach_id'][i][0]
                        node_id = dset['nodes']['node_id'][i][0]
            dset.close()
        except Exception as e:
            if not catch_exceptions:
                raise e
        dset = None
        if return_node:
            return node_id
        else:
            return reach_id
    
    def find_route_to_next_downstream_main_stem(self, reach_id : int, level : int) -> List[int]:
        """Internal: Returns a list of reaches on the way to the next main stem of the given Pfaffstetter level."""
        passed_reaches = []
        pfaff_2_code = str(reach_id)[:2]
        nc_path = self.get_dset_path(pfaff_2_code=pfaff_2_code)
        dset = self.get_dset(nc_path)

        while True:
            reach_index = np.where(dset['reaches']['reach_id'][:]==reach_id)[0][0]
            passed_reaches.append(reach_id)
            if int(str(reach_id)[level-1]) % 2 != 0:
                break
            if dset['reaches']['n_rch_down'][reach_index] == 0:
                raise SWORDConnectionException('End of Stream reached without passing query_reach')
            reach_id = dset['reaches']['rch_id_dn'][0][reach_index]
            if reach_id in passed_reaches:
                # Break loop in rare cases. e.g Wetlands
                break
        
        dset.close()
        dset = None
        return passed_reaches

    def find_route_along_main_stem(self, other_reach_id : int, reach_id : Optional[int] = None) -> List[int]:
        r"""Internal: Returns all passed reaches between the given reach_ids along the main stem. Function with overhead, but provides better results than find_route_along_main_stem_old"""
        reach_id = self.reach_id if (self.reach_id is not None and reach_id is None) else reach_id
        if reach_id is None:
            raise SWORDException('Reach_id argument or instance parameter required.')
        b1, b2 = str(other_reach_id)[0:6], str(reach_id)[0:6]
        if b1[:2] != b2[:2]:
            raise SWORDException(f'Reaches {other_reach_id} and {reach_id} not in same region.')
        pfaff_2_code = str(reach_id)[:2]
        nc_path = self.get_dset_path(pfaff_2_code=pfaff_2_code)
        dset = self.get_dset(nc_path)

        connected = False
        for fields in [{'n_field': 'n_rch_down', 'id_field': 'rch_id_dn' },{'n_field': 'n_rch_up', 'id_field': 'rch_id_up' }]:
            for start_id, end_id in [[reach_id,other_reach_id],[other_reach_id,reach_id]]:
                _reach_id = start_id
                passed_reaches = []
                while True:
                    reach_index = np.where(dset['reaches']['reach_id'][:]==_reach_id)[0][0]
                    passed_reaches.append(_reach_id)
                    if _reach_id == end_id:
                        connected = True
                        break
                    if dset['reaches'][fields['n_field']][reach_index] == 0:
                        break
                    _reach_id = dset['reaches'][fields['id_field']][0][reach_index]
                    if _reach_id in passed_reaches:
                        # Break loop in rare cases. e.g Wetlands
                        break
                if connected:
                    break
            if connected:
                break
        logger.debug(passed_reaches)
        dset.close()
        dset = None
        if connected:
            return passed_reaches
        else:
            raise SWORDConnectionException(f'No Route found between reach {reach_id} and {other_reach_id}')

    def get_all_downstream_reaches(self, reach_id : Optional[int] = None, recursive_passed : List[int] = []) -> List[int]:
        r"""Internal : Returns all downstream reach_ids"""
        passed_reaches = []
        reach_id = self.reach_id if (self.reach_id is not None and reach_id is None) else reach_id
        if reach_id is None:
            raise SWORDException('Reach_id argument or instance parameter required.')
        pfaff_2_code = str(reach_id)[:2]
        nc_path = self.get_dset_path(pfaff_2_code=pfaff_2_code)
        dset = self.get_dset(nc_path)

        df = pd.DataFrame({
            'reach_id' : dset['reaches']['reach_id'][:],
            'n_rch_down' : dset['reaches']['n_rch_down'][:],
            'n_rch_up' : dset['reaches']['n_rch_up'][:],
            'rch_id_dn_1' : dset['reaches']['rch_id_dn'][0][:],
            'rch_id_dn_2' : dset['reaches']['rch_id_dn'][1][:],
            'rch_id_dn_3' : dset['reaches']['rch_id_dn'][2][:],
            'rch_id_dn_4' : dset['reaches']['rch_id_dn'][3][:],
        })
        df.rch_id_dn_1.replace(0,np.nan,inplace=True)
        df.rch_id_dn_2.replace(0,np.nan,inplace=True)
        df.rch_id_dn_3.replace(0,np.nan,inplace=True)
        df.rch_id_dn_4.replace(0,np.nan,inplace=True)
        df.set_index('reach_id',drop=True,inplace=True)
        last_reach_id = reach_id
        dset.close()
        dset = None
        while True:
            n_rch_down = int(df.n_rch_down.loc[last_reach_id])
            if n_rch_down == 0:
                break
            elif n_rch_down > 1:
                reach_id = int(df.loc[reach_id,["rch_id_dn_1","rch_id_dn_2","rch_id_dn_3","rch_id_dn_4"]].min())
            else:
                reach_id = int(df.rch_id_dn_1.loc[last_reach_id])
            if reach_id in passed_reaches:
                break
            passed_reaches.append(reach_id)
            last_reach_id = reach_id
        return list(set(passed_reaches))

    def get_distance(self, lon1, lat1, lon2, lat2, reach_id1 : Optional[int] = None, reach_id2 : Optional[int] = None, reaches_between : Optional[List[int]] = None, max_cl_dist=1):
        if reach_id1 is None:
            reach_id1 = self.search_reach_id(lon1,lat1,max_dist=max_cl_dist)
        if reach_id2 is None:
            reach_id2 = self.search_reach_id(lon2,lat2,max_dist=max_cl_dist)
        if reach_id1 is None or reach_id2 is None:
            return None, None, reach_id1, reach_id2
        if reaches_between is None:
            reaches_between = self.find_route_along_main_stem(reach_id=reach_id1, other_reach_id=reach_id2)
        
        distance = 0
        for reach_id in reaches_between:
            if reach_id not in [reach_id1,reach_id2]:
                distance += SWORD_Reach(reach_id,self.sword_base_root,self.sword_version).length

        if reach_id1 > reach_id2:
            up_id, up_lon, up_lat, down_id, down_lon, down_lat = reach_id1, lon1, lat1, reach_id2, lon2, lat2
        else:
            up_id, up_lon, up_lat, down_id, down_lon, down_lat = reach_id2, lon2, lat2, reach_id1, lon1, lat1
        up_reach = SWORD_Reach(up_id,self.sword_base_root,self.sword_version)
        down_reach = SWORD_Reach(down_id,self.sword_base_root,self.sword_version)
        
        distance += up_reach.get_reach_position_of_closest_cl_point(up_lon,up_lat)
        distance += down_reach.length - down_reach.get_reach_position_of_closest_cl_point(down_lon,down_lat)
        return distance, reaches_between, reach_id1, reach_id2

    def check_connected(self, other_reach_id : int, reach_id : Optional[int] = None, max_tributaries : int = np.inf, brute_force : bool = False, silent : bool = False) -> bool:
        r"""
        
        Returns True if queried reach_id is connected with other_reach_id.
        Returns False if the reaches are not within the same river system or a dam/waterfall is between them.
        
        Parameters
        ----------
        other_reach_id : int
            reach_id of the reach to test the connection with.
        reach_id : int (Optional)
            reach_id of base reach for connection test. If the function is called from a SWORD_Reach instance and reach_id is not given, the instance reach_id is used.
        max_tributaries: int (Default: np.inf)
            If one of both reaches is not along the main stem of the river, the function returns False if the number of reaches between both reaches within the tributary is above the max_tributaries value.
        brute_force : bool (Default: False)
            if True, max_tributaries is ignored and a simple test is done if both reaches share a set of equal downstream reaches.
        silent : bool (Default: False)
            if True there will be no print or logging info output
        """
        reach_id = self.reach_id if (self.reach_id is not None and reach_id is None) else reach_id
        if reach_id is None:
            raise SWORDException('reach_id argument or instance parameter required.')
        reaches_between = []
        # b1, b2 = str(query_reach_id[0:6]), str(self.reach_id[0:6])
        # upstream_reach, downstream_reach = [query_reach_id, self.reach_id] if b1 > b2 else [self.reach_id, query_reach_id]
        logger.debug(f'Checking connection of reach {reach_id} and {other_reach_id}')
        try:
            if str(other_reach_id)[0:2] != str(reach_id)[0:2]:
                logger.debug('Not in same Level 2 basin')
                return False
            elif str(other_reach_id) == str(reach_id):
                # Same Reach
                logger.debug('Same Reach')
                return True
            elif str(other_reach_id)[:6] == str(reach_id)[:6]:
                # Same Level 6 Basin. Check for flow discontinuties
                logger.debug('Same Level 6 Basin')
                reaches_between = self.find_route_along_main_stem(other_reach_id, reach_id)
            elif brute_force:
                logger.warning('Brute Force Connection Check: Tributary constraints ignored.')
                downstream_reaches = self.get_all_downstream_reaches(reach_id)
                other_downstream_reaches = self.get_all_downstream_reaches(other_reach_id)
                shared = set(downstream_reaches) & set(other_downstream_reaches)
                if len(shared) == 0:
                    logger.debug('No shared reaches')
                    return False
                reaches_between = list(set(downstream_reaches).symmetric_difference(set(other_downstream_reaches)))
                
            elif str(other_reach_id)[0:2] == str(reach_id)[0:2]:
                # Same Level 2 basin
                logger.debug('Same Level 2 basin')
                n_tributaries = 0
                for level in [3,4,5,6]:
                    if (other_level := str(other_reach_id)[level-1]) == (this_level := str(reach_id)[level-1]):
                        logger.debug(f'Same Level {level} basin')
                        continue
                    else:
                        if int(other_level) % 2 != 0 and int(this_level) % 2 != 0:
                            # Both in main stream
                            logger.debug(f'Both in main stem')
                            reaches_between = self.find_route_along_main_stem(other_reach_id, reach_id)
                        else:
                            # Tributary
                            if int(other_level) % 2 == 0:
                                logger.debug(f'Other reach in tributary')
                                tributary_reaches = self.find_route_to_next_downstream_main_stem(other_reach_id, level)
                                n_tributaries += len(tributary_reaches) - 2
                                b1_estuary = tributary_reaches[-1]
                                reaches_between += tributary_reaches[1:-1]
                            else:
                                b1_estuary = other_reach_id
                            if int(this_level) % 2 == 0:
                                logger.debug(f'Reach in tributary')
                                tributary_reaches = self.find_route_to_next_downstream_main_stem(reach_id, level)
                                n_tributaries += len(tributary_reaches) - 2
                                b2_estuary = tributary_reaches[-1]
                                reaches_between += tributary_reaches[1:-1]
                            else:
                                b2_estuary = reach_id
                            if n_tributaries > max_tributaries:
                                logger.debug('Too many tributaries')
                                return False
                            reaches_between += self.find_route_along_main_stem(b1_estuary, b2_estuary)
                        break
            else:
                return False
            for _rid in reaches_between:
                if str(_rid)[-1] == '4' and _rid not in Utilities.waterfall_override[self.sword_version]:
                    if not silent:
                        logger.info('Waterfall or dam')
                    return False
            return True
        except SWORDConnectionException as e:
            if not silent:
                logger.critical(f'Reaches {reach_id} and {other_reach_id} not connected due to {str(e)}')
            return False

    def __getstate__(self):
        """Internal: hdf5 dset can not be pickled so it has to be removed prior pickling the reach object"""
        self.close_dset()
        return self.__dict__

    def close_dset(self) -> None:
        """Internal: closes the sword netcdf dataset """
        logger.debug('Closing SWORD h5 dset.')
        if self.dset is not None:
            self.dset.close()
            self.dset = None

    def get_dset(self, nc_path : Optional[str] = None):
        """Internal: Opens the sword netcdf dataset """
        nc_path = self.nc_path if (self.nc_path is not None and nc_path is None) else nc_path
        if nc_path is None:
            raise SWORDException('Netcdf Path argument or instance parameter required to load dset.')
        logger.debug(f'Loading {nc_path.name} as SWORD dset.')
        return h5py.File(str(nc_path), 'r')

    def get_dset_path(self, pfaff_2_code : Optional[int] = None, continent : Optional[str] = None, lon : Optional[float] = None, lat : Optional[float] = None) -> pathlib.Path:
        """ Internal: Returns the path of the SWORD dset located by the arguments. """
        continent = self.continent if (self.continent is not None and continent is None) else continent
        pfaff_2_code = self.pfaff_2_code if (self.pfaff_2_code is not None and pfaff_2_code is None) else pfaff_2_code
        if continent is None:
            if pfaff_2_code is None:
                if lon is None or lat is None:
                    raise SWORDException('Lon and Lat or Continent or Pfaffstetter Level 2 code argument resp. instance parameter required to find dset.')
                pfaff_2_code = self.get_pfaff2_code(lon,lat)
            continent = self.search_continent(pfaff_2_code)
        continent = continent.lower()
        nc_root = self.sword_root / 'netcdf'
        for nc_fn in [f.stem for f in nc_root.iterdir() if f.is_file()]:
            if f'{continent}_' in nc_fn:
                break
            else:
                nc_fn = None
        if nc_fn is None:
            raise SWORDException(f'SWORD File not found for continent {continent} in dir {str(nc_root)}')
        return nc_root / f'{nc_fn}.nc'

    def get_river_position_of_closest_cl_point(self, lon: float, lat : float, simple :bool = False, max_dist : float = np.inf) -> float:
        r"""
        Returns the along-stream distance of the closest centerline point to the river estuary.

        Parameters
        ----------
        lon : float
        
        lat : float
        
        simple : bool (Defult: False)
            If True, the dist_out parameter of the closest reach is used in combination with the node position within the reach.
            If False, the sum of the length of the passed reaches to the estuary are used in combination with the node position within the reach.
        
        max_dist : float (Default np.inf)
            maximum distance to centerline in km, if further returns None
        
        Returns
        ----------
        float
            Distance to river estuary [m]
        """
        reach_id = self.search_reach_id(lon=lon, lat=lat, max_dist=max_dist)
        if reach_id is None:
            return None
        reach = SWORD_Reach(reach_id,self.sword_base_root,self.sword_version)
        if simple:
            distance = reach.dist_out
        else:
            downstream_reaches = reach.get_all_downstream_reaches()
            distance = 0
            for _reach_id in downstream_reaches:
                if _reach_id != reach_id:
                    distance += SWORD_Reach(_reach_id,self.sword_base_root,self.sword_version).length
        distance += reach.get_reach_position_of_closest_cl_point(poi_lon=lon,poi_lat=lat)
        return distance

    def get_node_normal(self, node_id : int, pointing : str = 'right') -> LineString:
        r"""
        Returns the normal at the given node w.r.t. the up- and downstream nodes.

        Parameters
        ----------
        node_id : int
            ID of the query node
        pointing : str, default 'right'
            Choose from 'left', 'straight', 'right' to choose the pointing direction of the normal.  
            If 'left' or 'right' the base of the vector will be on the oppoiste site and the vector will have a length of 4 times the SWORD width.
            If 'straight' the base of the vector will be at the node and the lenght will be 100m.
        Returns
        ----------
        shapely.geometry.LineString
            Pointing to the right in downstream direction.
        """
        
        if pointing not in ['left', 'straight', 'right']:
            raise SWORDException(f"Pointing mode {pointing} unknown. Choose from 'left', 'straight', 'right'")
        node = SWORD_Node(node_id=node_id,sword_root=self.sword_base_root,sword_version=self.sword_version)
        return node.get_normal(pointing = pointing)    
    
    def export_to_shp(self, ids : List[int], path : pathlib.Path):
        schema = {
            'geometry': 'LineString',
            'properties': {'id': 'int'},
        }
        path = path.resolve()
        if not path.parent.exists():
            path.parent.mkdir()
        with fiona.open(path, mode='w', driver='ESRI Shapefile', schema=schema, crs=fiona.crs.from_epsg(4326)) as c:
            for rid in tqdm(ids, leave=False, desc='Exporting Reach geometries.'):
                reach = SWORD_Reach(rid)
                geom = LineString(reach.coords)
                c.write({
                    'geometry': mapping(geom),
                    'properties': {'id': int(rid)},
                })

class SWORD_Basin(SWORD_IF):
    """:meta private:"""
    pass

class SWORD_Node(SWORD_IF):
    def __init__(self,node_id : int, sword_root=None, sword_version=None):
        SWORD_IF.__init__(self,sword_root=sword_root,sword_version=sword_version)
        self.node_id = node_id
        self.pfaff_2_code = str(self.node_id)[:2]
        self.continent = self.search_continent(self.pfaff_2_code)

        self.nc_path = self.get_dset_path()
        self.nc_fn = self.nc_path.name
        self.dset = self.get_dset()

        self.node_index = np.where(self.dset['nodes']['node_id'][:]==self.node_id)[0][0]
        self.reach_id = self.dset['nodes']['reach_id'][self.node_index]
        self.reach_index = np.where(self.dset['reaches']['reach_id'][:]==self.reach_id)[0][0]
        self.slope = self.dset['reaches']['slope'][self.reach_index]/1000.
        self.lon, self.lat = self.dset['nodes']['x'][self.node_index], self.dset['nodes']['y'][self.node_index]
        self.width = self.dset['nodes']['width'][self.node_index]
        self.close_dset()
    
    def get_normal(self, pointing : str = 'right') -> LineString:
        if self.dset is None:
            self.dset = self.get_dset()
        assert pointing in ['left', 'straight', 'right'], f"Pointing mode {pointing} unknown. Choose from 'left', 'straight', 'right'"
        reach = SWORD_Reach(self.reach_id,self.sword_base_root,self.sword_version)
        short_indices = np.array(self.dset['nodes']['node_id'][:]) // 10
        this_node_position = int(str(self.node_id)[-3:-1])
        next_node_in_other_reach = False
        if this_node_position > 1:
            ds_node_index = np.where(short_indices==self.node_id // 10 - 1)
            if len(self.dset['nodes']['node_id'][ds_node_index]) > 0:
                ds_node_id = self.dset['nodes']['node_id'][ds_node_index][0]
            else:
                next_node_in_other_reach = True
        else:
            next_node_in_other_reach = True
        if next_node_in_other_reach:
            ds_reach_ids = reach.get_downstream_reach_ids()
            if len(ds_reach_ids) == 0:
                logger.info(f'No Downstream Reach for node : {self.node_id}. Tyring Substitution with node itself.')
                ds_node_id = None
            else:
                ds_reach = SWORD_Reach(ds_reach_ids[0],self.sword_base_root,self.sword_version)
                ds_short_id = int(f'{str(ds_reach.reach_id)[:-1]}{str(ds_reach.n_nodes).zfill(3)}')
                ds_node_index = np.where(short_indices==ds_short_id)
                if len(self.dset['nodes']['node_id'][ds_node_index]) > 0:
                    ds_node_id = self.dset['nodes']['node_id'][ds_node_index][0]
                else:
                    logger.info(f'No Downstream Node for node : {self.node_id}. Tyring Substitution with node itself.')
                    ds_node_id = None

        next_node_in_other_reach = False
        if this_node_position != reach.n_nodes:
            us_node_index = np.where(short_indices==self.node_id // 10 + 1)
            if len(self.dset['nodes']['node_id'][us_node_index]) > 0:
                us_node_id = self.dset['nodes']['node_id'][us_node_index][0]
            else:
                next_node_in_other_reach = True
        else:
            next_node_in_other_reach = True

        if next_node_in_other_reach:
            us_reach_ids = reach.get_upstream_reach_ids()
            if len(us_reach_ids) == 0:
                logger.info(f'No Upstream Reach for node : {self.node_id}. Tyring Substitution with node itself.')
                us_node_id = None
            else:
                us_reach = SWORD_Reach(us_reach_ids[0],self.sword_base_root,self.sword_version)
                us_short_id = int(f'{str(us_reach.reach_id)[:-1]}001')
                us_node_index = np.where(short_indices==us_short_id)
                if len(self.dset['nodes']['node_id'][us_node_index]) > 0:
                    us_node_id = self.dset['nodes']['node_id'][us_node_index][0]
                else:
                    logger.info(f'No Upstream Node for node : {self.node_id}. Tyring Substitution with node itself.')
                    us_node_id = None

        logger.debug(f'Upstream Node:   {us_node_id}')
        logger.debug(f'Downstream Node: {ds_node_id}')

        if us_node_id is None and ds_node_id is None:
            raise SWORDException(f'Can not calculate normal at node_id: {self.node_id}')
        elif us_node_id is None:
            us_node_id = self.node_id
            us_node_index = self.node_index
        elif ds_node_id is None:
            ds_node_id = self.node_id
            ds_node_index = self.node_index

        us_lon, us_lat = self.dset['nodes']['x'][us_node_index], self.dset['nodes']['y'][us_node_index]
        ds_lon, ds_lat = self.dset['nodes']['x'][ds_node_index], self.dset['nodes']['y'][ds_node_index]
        utm_zone = Conversions.convert_wgs_to_utm_zone(self.lon,self.lat)
        us_utm = Conversions.convert_wgs_geometry_to_utm(Point(us_lon,us_lat), utm_zone)
        ds_utm = Conversions.convert_wgs_geometry_to_utm(Point(ds_lon,ds_lat), utm_zone)
        this_utm = Conversions.convert_wgs_geometry_to_utm(Point(self.lon,self.lat), utm_zone)

        # Vector between up- and downstream nodes
        x, y = ds_utm.x - us_utm.x, ds_utm.y - us_utm.y
        norm = np.sqrt((x**2) + (y**2))
        x, y  = x/norm, y/norm
        if pointing == 'straight':
            normal = Conversions.convert_utm_geometry_to_wgs(LineString([this_utm, Point(this_utm.x + 100*x, this_utm.y + 100*y)]), utm_zone)
        else:
            # Orthogonal vector
            u, v = -y, x
            pointA = Point(this_utm.x + u*2*self.width, this_utm.y + v*2*self.width)
            
            u, v = y, -x
            pointB = Point(this_utm.x + u*2*self.width, this_utm.y + v*2*self.width)

            if pointing == 'right':
                normal = Conversions.convert_utm_geometry_to_wgs(LineString([pointA, pointB]), utm_zone)
            elif pointing == 'left':
                normal = Conversions.convert_utm_geometry_to_wgs(LineString([pointB, pointA]), utm_zone)
        
        self.close_dset()
        return normal 

class SWORD_Reach(SWORD_IF):
    @classmethod
    def get_pickle_path_from_reach_id(cls, reach_id : int, root : pathlib.Path) -> pathlib.Path:
        """
        classmethod to get the path from the given root for the given reach id.

        Parameters
        ---------
        reach_id : int
        root : pathlib.Path

        Returns
        ---------
        pathib.Path
        """
        path = root
        reach_id_str = str(reach_id)
        for x in reach_id_str[:6]:
            path = path / x
        return path / (reach_id_str + '.pickle')
    
    @classmethod
    def load(cls, reach_id : int, root : pathlib.Path):
        path = cls.get_pickle_path_from_reach_id(reach_id, root)
        with open(path, 'rb') as f:
            return pickle.load(f)

    def __init__(self,reach_id,sword_root=None,sword_version=None):
        SWORD_IF.__init__(self,sword_root=sword_root,sword_version=sword_version)
        assert True if isinstance(reach_id,(int, np.integer)) else float(reach_id).is_integer(), 'Reach ID must be an int like!'
        assert len(str(int(reach_id))) == 11, 'Reach ID must be 11 digits long!'
        self.reach_id = int(reach_id)
        self.pfaff_2_code = str(self.reach_id)[:2]
        self.continent = self.search_continent(self.pfaff_2_code)
        self.shp_path = self.get_shp_path()
        self.nc_path = self.get_dset_path()
        self.nc_fn = self.nc_path.name
        self.dset = self.get_dset()

        # Load SWORD Data
        self.reach_index = np.where(self.dset['reaches']['reach_id'][:]==self.reach_id)[0][0]

        if 'river_name' in self.dset['reaches'].keys():
            self.river_name = self.dset['reaches']['river_name'][self.reach_index]
        else:
            self.river_name = None

        self.dist_out = self.dset['reaches']['dist_out'][self.reach_index]
        self.slope = self.dset['reaches']['slope'][self.reach_index]/1000.
        self.dist_out = self.dset['reaches']['dist_out'][self.reach_index]
        self.width = self.dset['reaches']['width'][self.reach_index]
        self.width_std = np.sqrt(self.dset['reaches']['width_var'][self.reach_index])
        self.n_nodes = self.dset['reaches']['n_nodes'][self.reach_index]
        self.wse = self.dset['reaches']['wse'][self.reach_index]
        
        xmax = self.dset['reaches']['x_max'][self.reach_index]
        ymax = self.dset['reaches']['y_max'][self.reach_index]
        xmin = self.dset['reaches']['x_min'][self.reach_index]
        ymin = self.dset['reaches']['y_min'][self.reach_index]
        self.simple_center = [np.mean([xmax,xmin]),np.mean([ymax,ymin])]
        self.bbox = [xmin,ymin,xmax,ymax]
        self.coords = self.get_coordinates_from_shapefile()
        self.length = self.get_reach_length()
        # self.reference_point = self.get_reference_coord()
        self.upstream_reach_ids = self.get_upstream_reach_ids()
        self.downstream_reach_ids = self.get_downstream_reach_ids()
        self.aoi = self.get_aoi(std_multiplicator=0)
        self.jrc_aoi = None

        ## Advanced
        self.intersecting_passes = None
        self.intersected_misions = []
        self.icesat2_data = None
        self.icesat2_slope_data = None
        self.extended_icesat2_aoi = False
        self.icesat2_slope_data_no_conf_th = None
        self.mean_along_slope = None
        self.mean_across_slope = None
        self.mean_mixed_slope = None
        self.features_by_mission = None
        self.close_dset()

    def get_pickle_path(self, root : pathlib.Path):
        path = self.get_pickle_path_from_reach_id(self.reach_id, root)
        return path

    def save(self, path : pathlib.Path):
        if not path.parent.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
        with open(path, 'wb') as f:
            pickle.dump(self,f)

    def get_shp_path(self, pfaff_2_code : Optional[int] = None, continent : Optional[str] = None):
        """ Internal: Returns the path of the SWORD shapefile containing the reach """
        continent = self.continent if (self.continent is not None and continent is None) else continent
        pfaff_2_code = self.pfaff_2_code if (self.pfaff_2_code is not None and pfaff_2_code is None) else pfaff_2_code
        if continent is None:
            if pfaff_2_code is None:
                pfaff_2_code = self.get_pfaff2_code(self.simple_center[0],self.simple_center[1])
            continent = self.search_continent(pfaff_2_code)
        shp_root = self.sword_root / 'shp' / continent.upper()
        for shp_fn in [f.stem for f in shp_root.iterdir() if f.is_file()]:
            if f'sword_reaches_hb{self.pfaff_2_code}' in shp_fn:
                break
            else:
                shp_fn = None
        if shp_fn is None:
            raise SWORDException(f'SWORD File not found for pfaff_2_code {self.pfaff_2_code} in dir {str(shp_root)}')
        return shp_root / f'{shp_fn}.shp'

    def get_reach_length(self):
        lon, lat = zip(*self.coords)
        distances = Utilities.spherical_distance(np.array(lon[1:]),np.array(lat[1:]),np.array(lon[:-1]),np.array(lat[:-1]))
        return np.sum(distances) * 1000

    def get_reach_position_of_closest_cl_point(self, poi_lon : float, poi_lat : float) -> float:
        """ get the distance from the reach outflow of the closests reach centerline point to the given coord """
        # coords are (mostly) ordered from downstream to upstream
        lon, lat = zip(*self.coords)
        distances_to_poi = Utilities.spherical_distance(np.array(lon),np.array(lat),poi_lon,poi_lat)
        distances = np.insert(Utilities.spherical_distance(np.array(lon[1:]),np.array(lat[1:]),np.array(lon[:-1]),np.array(lat[:-1])),0,0) * 1000.
        closest_coord_index = distances_to_poi.argmin()
        
        if self.get_reach_position_of_closest_node(*self.coords[0]) < self.get_reach_position_of_closest_node(*self.coords[-1]):
            # typical order from downstream to upstream
            position = np.sum(distances[:closest_coord_index])
        else:
            # rare opposite order
            position = np.sum(distances) - np.sum(distances[:closest_coord_index])
        return position

    def get_coordinates_from_shapefile(self):
        """ Returns list of centerline coordinates """
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataSource = driver.Open(str(self.shp_path), 0)
        layer = dataSource.GetLayer()
        layer.SetAttributeFilter(f"reach_id ILIKE '{self.reach_id}'")
        feature = None
        for feature in layer:
            if int(feature.GetField("reach_id")) == self.reach_id:
                break
        if feature == None:
            raise SWORDException(f'No centerline shapefile feature for reach {self.reach_id}')
        geom = feature.GetGeometryRef()
        coords = []
        for i in range(geom.GetPointCount()):
            g = geom.GetPoint(i)
            coords.append(g[0:2])
        driver = None
        dataSource = None
        layer = None
        return coords

    def get_downstream_reach_ids(self):
        """ Returns downstream reach ids given in netcdf dataset """
        if self.dset is None:
            self.dset = self.get_dset()
        reach_ids = []
        for i in range(self.dset['reaches']['n_rch_down'][self.reach_index]):
            r_id = self.dset['reaches']['rch_id_dn'][i][self.reach_index]
            if r_id != 0:
                reach_ids.append(r_id)
        self.close_dset()
        return reach_ids

    def get_upstream_reach_ids(self):
        """ Returns upstrea, reach ids given in netcdf dataset """
        if self.dset is None:
            self.dset = self.get_dset()
        reach_ids = []
        for i in range(self.dset['reaches']['n_rch_up'][self.reach_index]):
            r_id = self.dset['reaches']['rch_id_up'][i][self.reach_index]
            if r_id != 0:
                reach_ids.append(r_id)
        self.close_dset()
        return reach_ids

    def get_aoi(self,buffersize=None, std_multiplicator = 1):
        """ get a coarse aoi of the reach buffered by the given parameter or the width plus four times the width standard deviation"""
        if buffersize is None:
            buffersize = self.width + std_multiplicator * self.width_std
        
        centerline_wgs84 = LineString(self.coords[1:-1])
        epsg_utm = Conversions.convert_wgs_to_utm_zone(geometry=centerline_wgs84)
        centerline_utm = Conversions.convert_wgs_geometry_to_utm(centerline_wgs84,epsg_utm)

        centerline_utm = centerline_utm.simplify(tolerance=100)
        aoi_utm = centerline_utm.buffer(buffersize, cap_style=2)
        aoi_wgs84 = Conversions.convert_utm_geometry_to_wgs(aoi_utm, epsg_utm)
        if aoi_wgs84.geom_type != 'Polygon':
            aoi_wgs84 = aoi_wgs84.convex_hull
        return list(aoi_wgs84.exterior.coords)

    def get_node_coords(self):
        """Returns a list of coords of the SWORD nodes"""
        if self.dset is None:
            self.dset = self.get_dset()
        nodes_indexes = np.where(self.dset['nodes']['reach_id'][:]==self.reach_id)[0]
        node_coords = [[self.dset['nodes']['x'][i],self.dset['nodes']['y'][i]] for i in nodes_indexes]
        self.close_dset()
        return node_coords
    
    def get_reach_position_of_closest_node(self, poi_lon : float, poi_lat : float) -> float:
        """ get the distance from the reach outflow of the closests reach node to the given coord """
        # nodes are ordered from downstream to upstream
        if self.dset is None:
            self.dset = self.get_dset()
        nodes_indexes = np.where(self.dset['nodes']['reach_id'][:]==self.reach_id)[0]
        lon = [self.dset['nodes']['x'][i] for i in nodes_indexes]
        lat = [self.dset['nodes']['y'][i] for i in nodes_indexes]
        df = pd.DataFrame({'lon': lon, 'lat':lat}, index = nodes_indexes).sort_index()
        closest_node_index = df.apply(lambda x: Utilities.spherical_distance(poi_lon,poi_lat,x.lon,x.lat),axis=1).idxmin()
        df.loc[:,'last_lon'] = [np.nan] + list(lon[:-1])
        df.loc[:,'last_lat'] = [np.nan] + list(lat[:-1]) 
        # spherical_distance returns distance in km!
        df.loc[:,'distance'] = df.apply(lambda x: Utilities.spherical_distance(x.lon,x.lat,x.last_lon,x.last_lat) * 1000,axis=1)
        df.distance.iloc[0] = 0
        position = df.loc[0:closest_node_index,'distance'].sum()
        self.close_dset()
        return position

    def contains_coord(self, lon : float, lat : float) -> bool:
        if Polygon(self.aoi).contains(Point(lon,lat)):
            return True
        else:
            return False

    def get_crossing_passes(self, dahiti_instance, mva_instance, ee_instance = None, focus_mission_groups : Optional[List[str]] = None, skip_mission_groups : Optional[List[str]] = None):
        """
        Intersects reach AOI with MVA orbits and returns crossing missions/passes/cycles.
        Returns the intersecting passes considering the given filters focus_mission_groups and skip_mission_groups.
        Appends intersecting passes to self.intersecting_passes.

        Parameters
        ----------
        dahiti_instance : Instance of dgfi_if.DAHITI_IF Class

        mva_instance : Instance of dgfi_if.MVA_IF Class

        ee_instance : Instance of dgfi_if.EE_IF Class (optional)
            If given, the JRC Water Occurence > 80 will be used as AOI.
            Otherwise the standard SWORD AOI will be used
        
        focus_mission_groups : List of strings (optional)
            list of keys from the dgfi_if.Utilities.mission_groups dictionary to be processed. If not given all will be processed except the skip_mission_groups.  
            ['jason', 'envisat', 'saral', 'sentinel3a', 'sentinel3b', 'cryosat2', 'icesat', 'icesat2_gt1l', 'icesat2_gt1r', 'icesat2_gt2l', 'icesat2_gt2r', 'icesat2_gt3l', 'icesat2_gt3r']
        
        skip_mission_groups : List of strings (optional)
            list of keys from the Misc.mission_groups dictionary to be skipped

        Returns
        ----------
        Pandas DataFrame
            intersecting passes considering the given filters focus_mission_groups and skip_mission_groups
        """
        aoi = Polygon(self.aoi)
        min_lon, min_lat, max_lon, max_lat = aoi.bounds

        node_centerline = self.get_node_coords()
        if ee_instance is not None:
            if len(node_centerline) < 2:
                self.intersecting_passes = pd.DataFrame()
                raise SWORDException('The Centerline requires at least 2 points.')
            self.jrc_aoi = ee_instance.get_JRC_aoi(centerline=node_centerline, coarse_aoi=self.aoi)
            aoi = self.jrc_aoi

        cryosat_modes = Utilities.get_cryosat_modes(*self.simple_center)

        if focus_mission_groups is None:
            if skip_mission_groups is None:
                skip_mission_groups = self.intersected_misions
            else:
                skip_mission_groups.append(self.intersected_misions)

        label_data = [(label, data) for label, data in Utilities.mission_groups.items() if all([focus_mission_groups is None or label in focus_mission_groups, skip_mission_groups is None or label not in skip_mission_groups])]
        intersecting_passes = []
        for mission_group_label, mission_group_data in tqdm(label_data, desc=f'Searching passes intersecting Reach {self.reach_id}', leave=False):
            passes_df = dahiti_instance.get_passes_in_bbox(mission=mission_group_data['passlocator'], min_lon=min_lon, min_lat=min_lat, max_lon=max_lon, max_lat=max_lat, convert_cryosat=False)
            if passes_df.empty:
                continue
            pass_filter = passes_df['pass_nr'].unique()
            for mva_mission in tqdm(mission_group_data['mva'], desc=f'Iterating {mission_group_label} MVA Missions', leave=False):
                if mva_mission == 'cryosatD_sar_hf' and 'SAR' not in cryosat_modes:
                    continue
                if mva_mission == 'cryosatD_sin_hf' and 'SIN' not in cryosat_modes:
                    continue
                pass_features = mva_instance.get_passes(mva_mission, pass_filter=pass_filter)
                
                for _feature in tqdm(pass_features, desc='Checking Intersection', leave=False):
                    _feature['geometry'] = wkt.loads(_feature['geometry'])
                    intersection = aoi.intersection(Utilities.validate_geometry(_feature['geometry']))
                    if intersection.is_empty:
                        continue
                    else:
                        _feature['intersection'] = intersection
                        intersecting_passes.append(_feature)
                        _feature = None
            if mission_group_label not in self.intersected_misions:
                self.intersected_misions.append(mission_group_label)
        
        intersecting_passes = pd.DataFrame(intersecting_passes)
        if self.intersecting_passes is None:
            self.intersecting_passes = intersecting_passes
        else:
            self.intersecting_passes = pd.concat([self.intersecting_passes,intersecting_passes])
            self.intersecting_passes = self.intersecting_passes.drop_duplicates(subset=['mission','pass_nr','cycle_nr']).reset_index(drop=True)
        
        return intersecting_passes

    def get_icesat2_data(self, mva_instance = None, ee_instance = None, extend_aoi : bool = False, debug : Union[bool,str] = False, silent : bool = False):
        """
        Parameters
        ----------
        mva_instance : Instance of dgfi_if.MVA_IF Class
            Required in case the passes have to be intersected and extracted from mva.

        ee_instance : Instance of dgfi_if.EE_IF Class (optional)
            In case the passes have to be intersected and extracted from mva:
            If given, the JRC Water Occurence > 80 will be used as AOI.
            Otherwise the standard SWORD AOI will be used

        extend_aoi : bool (Default: False)
            If true extends aoi by one reach downstream and updstream

        debug : bool
        """
        if mva_instance is None:
            from .MVA_IF import MVA_IF
            mva_instance = MVA_IF()

        icesat2_data = []
        node_centerline = self.get_node_coords()
        polygon = Polygon(self.aoi)
        if self.features_by_mission is None:
            if extend_aoi:
                if len(self.upstream_reach_ids) == 1 and len(self.upstream_reach_ids) == 1:
                    if self.check_connected(self.upstream_reach_ids[0], self.downstream_reach_ids[0], silent=True):
                        upstream_reach = SWORD_Reach(self.upstream_reach_ids[0])
                        downstream_reach = SWORD_Reach(self.downstream_reach_ids[0])
                        node_centerline = upstream_reach.get_node_coords() + self.get_node_coords() + downstream_reach.get_node_coords()
                        polygon = unary_union([Polygon(upstream_reach.aoi),Polygon(self.aoi),Polygon(upstream_reach.aoi)])
                        self.extended_icesat2_aoi = True
                        if polygon.geom_type == 'MultiPolygon':
                            polygon = polygon.convex_hull
                    else:
                        logger.critical(f'AOI of reach {self.reach_id} can not be extended.')
                else:
                    logger.critical(f'AOI of reach {self.reach_id} can not be extended.')
            else:
                self.extended_icesat2_aoi = False

            if ee_instance is not None:
                if len(node_centerline) < 2:
                    self.intersecting_passes = pd.DataFrame()
                    raise SWORDException('The Centerline requires at least 2 points.')
                ee_polygon = ee_instance.get_JRC_aoi(centerline=node_centerline, coarse_aoi=polygon, occ_th=50)
                if not ee_polygon.is_empty:
                    polygon = ee_polygon
            try:
                if polygon.geom_type == 'MultiPolygon':
                    polygon = polygon.convex_hull
                features_by_mission = mva_instance.get_icesat2_data_by_polygon(polygon,silent=silent)
            except MVAException as e:
                logger.critical(repr(e))
                logger.critical(f'No ICESat-2 Data for Reach {self.reach_id} due to MVA Error.')
                self.atl13_pass_data = pd.DataFrame()
                features_by_mission = {}
            self.features_by_mission = features_by_mission

        def rec_append(input):
            # Helper Function for nested data
            if isinstance(input,list):
                for item in input:
                    rec_append(item)
            else:
                icesat2_data.append(input)

        for features in tqdm(list(self.features_by_mission.values()), desc=f'Estimating ICESat-2 Heights and Along Slope for Reach {self.reach_id} by beam', leave=False):
            logger.debug(f'{len(features)} features')
            for feature in features:
                if feature['data'].empty:
                    logger.debug(f'empty feature')
                    continue
                data = self.icesat2_pass_helper(feature, debug = debug)
                if data is None:
                    logger.debug(f'noData feature')
                    continue
                if isinstance(data, list):
                    rec_append(data)
                else:
                    icesat2_data.append(data)
        self.icesat2_data = pd.DataFrame(icesat2_data)

    def icesat2_pass_helper(self, feature : Dict, max_median_deviation : float = 0.05, min_size_to_apply_window : int = 20, window_size : int = 7, max_slope : float = 0.0003, debug : Union[bool,str]= False):
        """Helper Function to reject outliers and convert along track slope to along centerline slope"""
        pass_geom = wkt.loads(feature['geometry'])
        centerline = LineString(self.coords)

        #######################################################################################
        ### First Intersection for splitting features if necessary
        intersection = pass_geom.intersection(centerline)
        if intersection.is_empty:
            nearest = nearest_points(pass_geom, centerline)
            ref_point_lon = nearest[0].x
            ref_point_lat = nearest[0].y
        elif intersection.geom_type == 'Point':
            ref_point_lon = intersection.x
            ref_point_lat = intersection.y
        elif intersection.geom_type == 'MultiPoint':
            distances = {}
            for i, point in enumerate(intersection.geoms):
                distances[i] = feature['data'].apply(lambda x: spherical_distance(point.x,point.y,x.lon,x.lat),axis=1)
            df = pd.DataFrame(distances)
            nearest_point = df.idxmin(axis=1)
            data = []
            for i in distances.keys():
                this_index = nearest_point[nearest_point == i]
                new_feature = deepcopy(feature)
                new_feature['data'] = new_feature['data'][new_feature['data'].index.isin(this_index.index)]
                if len(new_feature['data'].lon) < 2:
                    continue
                try:
                    new_feature['geometry'] = Utilities.linestring_helper(new_feature['data'].lon, new_feature['data'].lat)
                    _data = self.icesat2_pass_helper(feature=new_feature,max_median_deviation=max_median_deviation,min_size_to_apply_window=min_size_to_apply_window,window_size=window_size,debug=debug)
                except Utilities.MiscException as e:
                    _data = None
                except SWORDException as e:
                    _data = None
                if _data is not None:
                    data.append(_data)
            return data
        #######################################################################################

        atl13_pass_data = feature['data']
        epsg_utm = Utilities.convert_wgs_to_utm_zone(lon=ref_point_lon, lat=ref_point_lat)
        ydata = atl13_pass_data.elev - atl13_pass_data.geoh07
        xdata = atl13_pass_data.alongtrack_distance

        #######################################################################################
        ### Outlier Detection

        if ydata.shape[0] > min_size_to_apply_window:
            window_median = ydata.rolling(window_size,min_periods=int(np.floor(window_size/2)),center=True).median()
            median_flags = (window_median - ydata).abs() <= max_median_deviation # False : AMD Outlier
        else:
            median_flags = (ydata - ydata.median()) <= max_median_deviation # False : AMD Outlier
        
        type_flags = True #atl13_pass_data.water_body_type.isin([2,5,6]) # Allowed Values: Reservoir, River, Estuary -> False : Type Outlier
        cloud_flags = atl13_pass_data.cloud_flag_asr_atl09.isin([0,1,2,3]) # -> False: Cloud Outlier
        ice_flags = atl13_pass_data.snow_ice_atl09.isin([0,1]) # -> False: Ice Outlier
        flags = type_flags & median_flags & cloud_flags & ice_flags # False : Type, Cloud, Ice, and AMD Outlier
        
        atl13_pass_data.loc[:,"type_flag"] = type_flags
        atl13_pass_data.loc[:,"median_flag"] = median_flags
        atl13_pass_data.loc[:,"cloud_flag"] = cloud_flags
        atl13_pass_data.loc[:,"ice_flags"] = ice_flags
        
        xdata_ice_outlier = xdata.loc[~ice_flags]
        ydata_ice_outlier = ydata.loc[~ice_flags]

        amd_xdata_outlier = xdata.loc[~flags] # For Debug Plot
        amd_ydata_outlier = ydata.loc[~flags] # For Debug Plot
        xdata = xdata.loc[flags]
        ydata = ydata.loc[flags]
        
        if xdata.shape[0] < 2:
            logger.debug('No Data after AMD and Type Outlier removal')
            return None
        
        clusters = (xdata.diff() > 500).cumsum()
        cluster_flags = clusters == clusters.value_counts().idxmax() # False: Cluster Outlier (Only longest Cluster is used)

        cluster_xdata_outlier = xdata.loc[~cluster_flags] # For Debug Plot
        cluster_ydata_outlier = ydata.loc[~cluster_flags] # For Debug Plot
        xdata = xdata.loc[cluster_flags]
        ydata = ydata.loc[cluster_flags]

        atl13_pass_data.loc[:,"cluster_flag"] = atl13_pass_data.index.isin(xdata.index) #cluster_flags

        if xdata.shape[0] < 2:
            logger.debug('No Data after Cluster Outlier removal')
            return None

        #######################################################################################
        ###  Position and Angle Determination
        pass_geom = Utilities.validate_geometry(LineString(atl13_pass_data[atl13_pass_data.cluster_flag.fillna(False)].wgs_coord.to_list()))
        intersection = pass_geom.intersection(centerline)
        if intersection.is_empty:
            nearest = nearest_points(pass_geom, centerline)
            ref_point_lon = nearest[0].x
            ref_point_lat = nearest[0].y
        elif intersection.geom_type == 'Point':
            ref_point_lon = intersection.x
            ref_point_lat = intersection.y
        elif intersection.geom_type == 'MultiPoint':
            raise SWORDException('Feature should be split')
        atl13_pass_data.loc[:,'ref_distance'] = Utilities.great_circle_dist(atl13_pass_data.lon.to_numpy(), atl13_pass_data.lat.to_numpy(), ref_point_lon, ref_point_lat) * 1000
        x0 = atl13_pass_data.alongtrack_distance.loc[atl13_pass_data.ref_distance.idxmin()]
        postion = self.get_reach_position_of_closest_cl_point(ref_point_lon, ref_point_lat)
        node_id = self.search_reach_id(lon=ref_point_lon, lat=ref_point_lat, return_node=True)
        normal = self.get_node_normal(node_id, pointing='straight')
        normal_utm = Utilities.convert_wgs_geometry_to_utm(normal, epsg_utm)
        normal_coords = np.asarray(Utilities.get_geom_coordinates(normal_utm))
        pass_geom_utm = Utilities.convert_wgs_geometry_to_utm(pass_geom, epsg_utm)
        pass_coords = np.asarray(Utilities.get_geom_coordinates(pass_geom_utm))
        
        normal_vector = normal_coords[-1] - normal_coords[0]
        pass_vector = pass_coords[-1] - pass_coords[0]
        u, v = pass_vector, normal_vector
        u_norm = np.sqrt(sum(u*u))
        v_norm = np.sqrt(sum(v*v))
        dot = np.dot(u, v)

        #######################################################################################
        ### This Block handles pass height estmate using SVR based on DAHITI approach and SVR outlier rejection
        limit = ((xdata.max() - xdata.min()) * max_slope) / 2
        limit = .05
        distances = np.abs(xdata - x0)
        distance_weights = 1./np.where(distances==0,1,distances)
        xscaler = StandardScaler() # SVR works better when data is scaled
        yscaler = StandardScaler()
        xscaler.fit(xdata.to_numpy().reshape(-1, 1))
        yscaler.fit(ydata.to_numpy().reshape(-1, 1))
        xscaled = xscaler.transform(xdata.to_numpy().reshape(-1, 1))
        yscaled = yscaler.transform(ydata.to_numpy().reshape(-1, 1))
        with warnings.catch_warnings():
            warnings.filterwarnings(action='ignore',message='Liblinear failed to converge, increase the number of iterations.')
            svr = LinearSVR(random_state=0,max_iter=1e6)
            svr.fit(xscaled,np.ravel(yscaled),sample_weight=distance_weights)
        y_predict = yscaler.inverse_transform(svr.predict(xscaled))
        y_deviation = (ydata - y_predict).abs()
        svr_flags = y_deviation <= limit
        y_predict = y_predict[svr_flags]
        if y_predict.size > 0:
            # svr_elev = np.nanmean(y_predict)
            # distances = distances[svr_flags]
            svr_elev = np.average(ydata[svr_flags],weights=distance_weights[svr_flags])
            median_elev = np.median(ydata)
        else:
            svr_elev = np.nan
            median_elev = np.nan

        svr_xdata_outlier = xdata.loc[~svr_flags] # For Debug Plot
        svr_ydata_outlier = ydata.loc[~svr_flags] # For Debug Plot
        xdata = xdata[svr_flags]
        ydata = ydata[svr_flags]
        width = np.max(xdata) - np.min(xdata)
        atl13_pass_data.loc[:,'svr_flag'] = atl13_pass_data.index.isin(xdata.index)
        atl13_pass_data.loc[:,'rejected'] = ~(atl13_pass_data.svr_flag & atl13_pass_data.cluster_flag & atl13_pass_data.median_flag & atl13_pass_data.type_flag)

        if xdata.shape[0] < 2:
            if debug:
                logger.debug('No Data after SVR Outlier removal')
            return None

        #######################################################################################
        ### This Block handles the conversion from atl13 along track slope to along river centerline slope
        dot_crossing_angle_deg = np.rad2deg(np.arccos(dot/(u_norm * v_norm)))
        if dot != 0 and xdata.size > 1: # orthogonal flying
            tinv = lambda p, df: abs(stats.t.ppf(p/2, df))
            ts = tinv(0.05, len(xdata)-2)
            slope, intercept, r, p, se = stats.linregress(xdata, ydata)
            slope_confidence = ts*se 
            dh = -1 * (u_norm * slope)
            proj_of_u_on_v = (dot/v_norm**2)*v
            l = np.sqrt(sum(proj_of_u_on_v**2))
            adjusted_custom_along_track_slope = (dh/l) * np.sign(dot)
        else:
            adjusted_custom_along_track_slope = np.nan
            slope, intercept, r, p, se = np.nan,np.nan,np.nan,np.nan,np.nan
            slope_confidence = np.nan
        
        #######################################################################################

        date_object = julianDayDate(np.nanmean(atl13_pass_data.jday))
        date = datetime.datetime(year=date_object["year"], month=date_object["month"], day=date_object["day"], hour=date_object["hour"], minute=date_object["minute"])

        if debug is True or debug == date.strftime('%Y-%m-%d'):
            plt.rc('font', size=14)  
            fig, ax = plt.subplots(figsize=(20,10))
            plt.sca(ax)
            cloud_flag_asr_atl09 = atl13_pass_data[atl13_pass_data.svr_flag].cloud_flag_asr_atl09.to_numpy()
            sc = plt.scatter(xdata, ydata, c=cloud_flag_asr_atl09, label='Valid ATL13 Observations', cmap=matplotlib.colors.ListedColormap([tum.blue,tum.light_blue2,tum.light_blue,tum.light_grey,tum.dark_grey,tum.orange]), vmax=6)
            ylim = list(plt.ylim())
            xlim = list(plt.xlim())
            plt.plot(amd_xdata_outlier, amd_ydata_outlier, '.',color=tum.grey, label='Type, Cloud, or AMD Outlier')
            plt.plot(svr_xdata_outlier, svr_ydata_outlier, '+',color=tum.grey, label='SVR Outlier')
            plt.plot(cluster_xdata_outlier, cluster_ydata_outlier, 'x',color=tum.grey, label='Cluster Outlier')
            plt.plot(xdata_ice_outlier, ydata_ice_outlier, '.',color='red', label='ICE Outlier')
            # plt.colorbar(sc)
            plt.plot(xdata, intercept + slope*xdata, color=tum.orange, label='Slope Fit')
            # plt.plot(xdata,y_predict,'b--',label="SVR fit")
            if svr_elev is not None:
                plt.axhline(svr_elev, label='Reference Point/Elevation',color=tum.grey)
            else:
                plt.axhline(median_elev, label='Median Height', linestyle='--')
            plt.axvline(x0,color=tum.grey)
            plt.xlabel('Along Track Distance [m]')
            plt.ylabel('Elevation [m]')
            date_str = date.strftime("%Y-%m-%d")
            plt.title(f'Reach {self.reach_id} Date {date_str} Beam {feature["beam"]}\nCrossing Angle: {dot_crossing_angle_deg:.2f} Along Track Slope: {adjusted_custom_along_track_slope*1e6:.0f} +- {slope_confidence*1e6:.0f} mm/km')
            quantiles = (atl13_pass_data.elev - atl13_pass_data.geoh07).quantile([0.1,0.99]).to_list()
            if ylim[0] > quantiles[0]:
                ylim[0] = quantiles[0]
            if ylim[1] < quantiles[1]:
                ylim[1] = quantiles[1]
            plt.ylim(ylim)
            # plt.xlim(xlim)
            plt.legend(loc='upper left')
            plt.show()
            plt.close()
            _df = gpd.GeoDataFrame(atl13_pass_data,geometry=atl13_pass_data.wgs_coord.apply(Point))
            _df = _df.set_crs(epsg=4326)
            _df = _df.to_crs(epsg=3857)
            _centerline = gpd.GeoSeries(centerline)
            _centerline = _centerline.set_crs(epsg=4326)
            _centerline = _centerline.to_crs(epsg=3857)
            ax = _df[~_df.rejected].plot(figsize=(6, 6),color=tum.light_blue2, label="Valid ATL13 Observations")
            _centerline.plot(ax=ax,linewidth=2, label="centerline",color=tum.orange)

            _aoi = gpd.GeoSeries(Polygon(self.aoi))
            _aoi = _aoi.set_crs(epsg=4326)
            _aoi = _aoi.to_crs(epsg=3857)
            _aoi.plot(ax=ax,linewidth=2, label="centerline",color=tum.blue,alpha=0.5)

            # ax.axis('equal')
            ymin, ymax = plt.ylim()
            xmin, xmax = plt.xlim()
            xamp = xmax - xmin
            yamp = ymax - ymin
            buffer = np.abs((xamp - yamp) / 2.)
            if xamp > yamp:
                ymax += buffer
                ymin -= buffer
            else:
                xmax += buffer
                xmin -= buffer

            bounds = gpd.GeoSeries(box(xmin,ymin,xmax,ymax)).set_crs(3857)
            bounds.plot(ax=ax,facecolor="none",edgecolor="none")
            bounds = bounds.to_crs(4326).iloc[0]
            # plt.ylim(np.array(plt.ylim()) + np.array([-1e4,1e4]))
            # plt.xlim(np.array(plt.xlim()) + np.array([-1e4,1e4]))
            try:
                # cx.add_basemap(ax, crs=_df.crs, source=cx.providers.Esri.WorldImagery)
                from .import EE_IF
                ee = EE_IF()
                url = ee.get_s2_mosaic_url(date=date,poi=bounds,daysdelta=21)
                cx.add_basemap(ax, crs=_df.crs, source=url)
            except:
                ...
            ax.add_artist(ScaleBar(2,location="lower right"))
            plt.legend()
            plt.xticks([], [])
            plt.yticks([], [])
            plt.show()
            plt.close()
                
        with warnings.catch_warnings():
            warnings.filterwarnings(action='ignore')
            return {
                'date': date.date(),
                'timestamp': date,
                'mission': feature['mission'],
                'position': postion,
                'lon': ref_point_lon,
                'lat': ref_point_lat,
                'median_elev': median_elev,
                'svr_elev' : svr_elev,
                'elev' : svr_elev if svr_elev is not None else median_elev,
                'elev_std': atl13_pass_data.elev[~atl13_pass_data.rejected].std(),
                'atl13_slope': atl13_pass_data.atl13_slope[~atl13_pass_data.rejected].median(),
                'min_depth': atl13_pass_data.depth[~atl13_pass_data.rejected].min(),
                'max_depth': atl13_pass_data.depth[~atl13_pass_data.rejected].max(),
                'median_depth': atl13_pass_data.depth[~atl13_pass_data.rejected].median(),
                'water_body_ids' : atl13_pass_data.water_body_id[~atl13_pass_data.rejected].unique(),
                'water_body_types' : atl13_pass_data.water_body_type[~atl13_pass_data.rejected].unique(),
                'hf_data': atl13_pass_data,
                'custom_along_track_slope' : slope,
                'custom_along_track_slope_err': se,
                'custom_along_track_slope_abs_rvalue': np.abs(r),
                'custom_along_track_slope_pvalue': p,
                'adjusted_custom_along_track_slope' : adjusted_custom_along_track_slope,
                'dot_crossing_angle_deg' : dot_crossing_angle_deg,
                'node_id' : node_id,
                'width' : width,
                'slope_confidence': slope_confidence,
                }
        
    def get_icesat2_slope(self, mva_instance = None, ee_instance = None, min_distance = 1000., debug=False, extend_aoi=False, extend_aoi_on_empty=False, silent=False, along_angle_th : Optional[float] = None, along_conf_th : Optional[float] = None):
        """
        Convenience function to extract ICESat2 Data from MVA and calculate water surfac slope.
        Sets the icesat2_data, icesat2_slope_data, mean_mixed_slope, mean_across_slope, and mean_along_slope attributes.

        Parameters
        ----------
        mva_instance : Instance of dgfi_if.MVA_IF Class
            Required in case the passes have to be intersected and extracted from mva.

        ee_instance : Instance of dgfi_if.EE_IF Class, optional
            In case the passes have to be intersected and extracted from mva:  
            If given, the JRC Water Occurence > 80 will be used as AOI.  
            Otherwise the standard SWORD AOI will be used (default is None)

        min_distance : float
            Minimum Distance betwen across beams to calculate across slope (default is 1000.)

        debug : bool
            prints debug messages and plots (default is False)
        
        extend_aoi : bool
            Extends AOI to upstream and downstream reaches (default is False)

        extend_aoi_on_empty : bool
            Extends AOI to upstream and downstream reaches if there is no data for the instance reach (default is False)

        silent : bool
            Suppresses any output (default is False)

        along_angle_th : float, optional
            Maximum angle argument for the along track outlier rejection (Recommended is 65, default is None)

        along_conf_th : float, optional
            Maximum slope fit confidence interval argument for the along track outlier rejection (Recommended is 300, default is None)

        Returns
        -------
        None
        """
        if mva_instance is None:
            from dgfi_if import MVA_IF
            mva_instance = MVA_IF(use_cache=False)

        if self.icesat2_data is None:
            self.get_icesat2_data(mva_instance, ee_instance, debug=debug, extend_aoi = extend_aoi, silent=silent)
        if self.icesat2_data.empty and extend_aoi_on_empty and not extend_aoi:
            self.get_icesat2_data(mva_instance, ee_instance, debug=debug, extend_aoi = True, silent=silent)
        if self.icesat2_data.empty:
            self.icesat2_slope_data = pd.DataFrame()
        else:
            self.icesat2_slope_data_no_conf_th = self.icesat2_slope_helper(min_distance = min_distance)
            self.icesat2_slope_data = self.icesat2_slope_helper(min_distance = min_distance, along_angle_th = along_angle_th, along_conf_th = along_conf_th)
            count = self.icesat2_slope_data.count()
            if extend_aoi_on_empty and 'mean_across' in count and 'mean_along' in count and count.mean_across < 0.25 * count.mean_along:
                self.get_icesat2_data(mva_instance, ee_instance, debug=debug, extend_aoi = True, silent=silent)
                self.icesat2_slope_data = self.icesat2_slope_helper(min_distance = min_distance, along_angle_th = along_angle_th, along_conf_th = along_conf_th)

        if "mean_along" in self.icesat2_slope_data:
            try:
                self.mean_along_slope = self.icesat2_slope_data.mean_along.median()
            except Exception as e:
                logger.warning('Cannot get along slope for reach {self.reach_id}')
                logger.exception(e, exc_info=True)
                self.mean_along_slope = np.nan
        else:
            self.mean_along_slope = np.nan
        if "mean_across" in self.icesat2_slope_data:
            try:
                self.mean_across_slope = np.average(self.icesat2_slope_data.mean_across.dropna(),weights=1./self.icesat2_slope_data.mean_across_std.dropna())
            except ZeroDivisionError:
                self.mean_across_slope = np.nanmean(self.icesat2_slope_data.mean_across.dropna())
            except Exception:
                logger.warning('Cannot get across slope for reach {self.reach_id}')
                logger.exception(e, exc_info=True)
                self.mean_across_slope = np.nan
        else:
            self.mean_across_slope = np.nan
        if "mean_along" in self.icesat2_slope_data and "mean_across" in self.icesat2_slope_data:
            self.mean_mixed_slope = self.icesat2_slope_data.mean_across.dropna().combine_first(self.icesat2_slope_data.mean_along.dropna()).median()
        else:
            self.mean_mixed_slope = np.nanmean([self.mean_across_slope, self.mean_along_slope])

    def icesat2_slope_helper(self, min_distance : float, along_angle_th : Optional[float] = None, along_conf_th : Optional[float] = None):
        mean_across_slope_by_date = {}
        mean_across_std_by_date = {}
        across_slopes_by_date = {}
        across_distances_by_date = {}
        mean_along_slope_by_date = {}
        across_x_angles_by_date = {}
        across_y_angles_by_date = {}
        mean_along_angle_by_date = {}
        along_slopes_by_date = {}
        water_body_types_by_date = {}
        water_body_ids_by_date = {}
        date_groups = dict(list(self.icesat2_data.groupby('date')))
        for date, df in date_groups.items():
            
            if len(df.index) > 1:
                ########### Across Slope
                cols = ['mission','position','elev','median_elev','elev_std','dot_crossing_angle_deg']
                cross = df[cols].merge(df[cols], how='cross')
                ### absolute distance only to remove duplicates to too close observations:
                cross['distance'] = np.abs(cross.position_x - cross.position_y)
                cross = cross[cross['distance'] > min_distance]
                cross = cross.drop_duplicates(subset="distance").reset_index(drop=True)
                if not cross.empty:
                    ###
                    ### Actual distance
                    cross['distance'] = cross.position_x - cross.position_y
                    cross['delta_h'] = cross.elev_x - cross.elev_y
                    cross['slope'] = cross.delta_h / cross.distance
                    cross['weight'] = 1./(cross.elev_std_x + cross.elev_std_y)
                    cross = cross[cross['slope'] >= 0]
                    try:
                        mean_across_slope_by_date[date] = np.average(cross.slope,weights=cross.weight)
                        mean_across_std_by_date[date] = np.average((cross.elev_std_x + cross.elev_std_y),weights=cross.weight)
                    except ZeroDivisionError:
                        mean_across_slope_by_date[date] = np.nanmean(cross.slope)
                        mean_across_std_by_date[date] = np.nanmean((cross.elev_std_x + cross.elev_std_y))
                    across_slopes_by_date[date] = cross['slope']
                    across_distances_by_date[date] = cross['distance']
                    across_x_angles_by_date[date] = cross['dot_crossing_angle_deg_x']
                    across_y_angles_by_date[date] = cross['dot_crossing_angle_deg_y']
            if len(df.index) > 0:
                ########### Along Slope
                def is_unique(s):
                    a = s.values
                    return (a[0] == a).all()
                
                if along_conf_th is not None and along_angle_th is not None:
                    df.loc[:,'conf_threshold'] = df['dot_crossing_angle_deg'].apply(Utilities.icesat2_along_get_conf_threshold_by_angle, args=(along_conf_th, along_angle_th))
                    df = df[df.slope_confidence * 1e6 < df.conf_threshold]
                    if df.empty:
                        continue
                
                p_index = df.adjusted_custom_along_track_slope >= 0
                df.loc[:,'reduced_angle_deg'] = df.dot_crossing_angle_deg.apply(lambda x: 180-x if x > 90 else x)
                if np.sum(p_index) > 0:
                    along_slopes_by_date[date] = df.adjusted_custom_along_track_slope
                    water_body_types_by_date[date] = is_unique(df.water_body_types.apply(pd.Series).stack().reset_index(drop=True))
                    water_body_ids_by_date[date] = is_unique(df.water_body_ids.apply(pd.Series).stack().reset_index(drop=True))
                    try:
                        mean_along_slope_by_date[date] = np.average(df.adjusted_custom_along_track_slope[p_index], weights=1./df.reduced_angle_deg[p_index])
                        mean_along_angle_by_date[date] = np.average(df.dot_crossing_angle_deg[p_index], weights=1./df.reduced_angle_deg[p_index])
                    except ZeroDivisionError:
                        mean_along_slope_by_date[date] = np.nanmean(df.adjusted_custom_along_track_slope[p_index])
                        mean_along_angle_by_date[date] = np.nanmean(df.dot_crossing_angle_deg[p_index])
                
        return pd.DataFrame({
            "mean_across" : pd.Series(mean_across_slope_by_date),
            "mean_along" : pd.Series(mean_along_slope_by_date),
            "across" : pd.Series(across_slopes_by_date),
            "along" : pd.Series(along_slopes_by_date),
            "mean_across_std" : pd.Series(mean_across_std_by_date),
            "mean_along_angle" : pd.Series(mean_along_angle_by_date),
            "across_distances" : pd.Series(across_distances_by_date),
            "across_angle_x" : pd.Series(across_x_angles_by_date),
            "across_angle_y" : pd.Series(across_y_angles_by_date),
            "unique_water_body_types" : pd.Series(water_body_types_by_date),
            "unique_water_body_ids" : pd.Series(water_body_ids_by_date),
        })