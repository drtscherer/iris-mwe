import netCDF4
import time
import pathlib
from typing import Union, Optional, List, Dict
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import pickle
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map, thread_map
from datetime import date
from dgfi_if import SWORD_Reach, SWORD_IF
from shapely.geometry import Point
tqdm.pandas()

SWORD_DOI_LOOKUP = {
    'v1r1' : '10.5281/zenodo.4917236',
    'v2r0' : '10.5281/zenodo.5643392',
}

def to_df(
        out_version : str,
        SWORD_version : str,
        root_dir : Union[str,pathlib.Path] = ...,
        out_dir : Union[str,pathlib.Path] = ".",
        max_workers : int = 20
        ):
    if not isinstance(root_dir, pathlib.Path):
        root_dir = pathlib.Path(root_dir)
    if not isinstance(out_dir, pathlib.Path):
        out_dir = pathlib.Path(out_dir)
    # root_dir = root_dir / ('reaches_' + SWORD_version)
    assert root_dir.exists(), f'Root dir "{str(root_dir)}" not found!'
    assert out_dir.exists(), f'Out dir "{str(root_dir)}" not found!'

    reach_ids = []

    for _, _, files in os.walk(root_dir):
        for f in files:
            reach_ids.append(int(f[0:-7]))

    def helper(reach_id):
        obj = SWORD_Reach.load(reach_id,root_dir)
        
        if "mean_along" in obj.icesat2_slope_data and "mean_across" in obj.icesat2_slope_data:
            mixed_slope_series = obj.icesat2_slope_data.mean_across.dropna().combine_first(obj.icesat2_slope_data.mean_along.dropna())
        elif 'mean_along' in obj.icesat2_slope_data:
            mixed_slope_series = obj.icesat2_slope_data.mean_along
        elif 'mean_across' in obj.icesat2_slope_data:
            mixed_slope_series = obj.icesat2_slope_data.mean_across
        else:
            mixed_slope_series = None

        return {
            'reach_id' : reach_id,
            'lon' : obj.simple_center[0],
            'lat' : obj.simple_center[1],
            
            'across_flag' : False if np.isnan(obj.mean_across_slope) else True,
            'along_flag' : False if np.isnan(obj.mean_along_slope) else True,
            'combined_flag' : False if np.isnan(obj.mean_mixed_slope) else True,
            
            'avg_across_slope' : obj.mean_across_slope * 1e6,
            'avg_along_slope' : obj.mean_along_slope * 1e6,
            'avg_combined_slope' : obj.mean_mixed_slope * 1e6,
            
            'min_across_slope' : obj.icesat2_slope_data.mean_across.min() * 1e6 if 'mean_across' in obj.icesat2_slope_data else np.nan,
            'min_along_slope' : obj.icesat2_slope_data.mean_along.min() * 1e6 if 'mean_along' in obj.icesat2_slope_data else np.nan,
            'min_combined_slope' : mixed_slope_series.min() * 1e6 if mixed_slope_series is not None else np.nan,
            
            'max_across_slope' : obj.icesat2_slope_data.mean_across.max() * 1e6 if 'mean_across' in obj.icesat2_slope_data else np.nan,
            'max_along_slope' : obj.icesat2_slope_data.mean_along.max() * 1e6 if 'mean_along' in obj.icesat2_slope_data else np.nan,
            'max_combined_slope' : mixed_slope_series.max() * 1e6 if mixed_slope_series is not None else np.nan,

            'std_across_slope' : obj.icesat2_slope_data.mean_across.std() * 1e6 if 'mean_across' in obj.icesat2_slope_data else np.nan,
            'std_along_slope' : obj.icesat2_slope_data.mean_along.std() * 1e6 if 'mean_along' in obj.icesat2_slope_data else np.nan,
            'std_combined_slope' : mixed_slope_series.std() * 1e6 if mixed_slope_series is not None else np.nan,
        
            'n_across_slope' : obj.icesat2_slope_data.mean_across.shape[0] if 'mean_across' in obj.icesat2_slope_data else 0,
            'n_along_slope' : obj.icesat2_slope_data.mean_along.shape[0] if 'mean_along' in obj.icesat2_slope_data else 0,
            'n_combined_slope' : mixed_slope_series.shape[0] if mixed_slope_series is not None else 0,
            
            'min_date_across_slope' : obj.icesat2_slope_data.mean_across.dropna().index.min() if 'mean_across' in obj.icesat2_slope_data else np.nan,
            'min_date_along_slope' : obj.icesat2_slope_data.mean_along.dropna().index.min() if 'mean_along' in obj.icesat2_slope_data else np.nan,
            'min_date_combined_slope' : mixed_slope_series.index.min() if mixed_slope_series is not None else np.nan,

            'max_date_across_slope' : obj.icesat2_slope_data.mean_across.dropna().index.max() if 'mean_across' in obj.icesat2_slope_data else np.nan,
            'max_date_along_slope' : obj.icesat2_slope_data.mean_along.dropna().index.max() if 'mean_along' in obj.icesat2_slope_data else np.nan,
            'max_date_combined_slope' : mixed_slope_series.index.max() if mixed_slope_series is not None else np.nan,

            'n_features' : sum([sum([not x['data'].empty for x in f]) for f in obj.features_by_mission.values()]) if obj.features_by_mission is not None else np.nan,
            'n_icesat_2_data' : obj.icesat2_data.shape[0] if obj.icesat2_data is not None else np.nan,
            'n_ice_flags' : obj.icesat2_data.hf_data.apply(lambda x: (~x.snow_ice_atl09.isin([0,1])).any()).sum() if obj.icesat2_data is not None and "hf_data" in obj.icesat2_data.columns else np.nan

            }
            
    data = thread_map(helper,reach_ids,desc='Raeding Reaches',max_workers=max_workers)
    df = pd.DataFrame.from_records(data)
    df.to_pickle(str(out_dir / ('IRIS_pandas_' + out_version + '.pkl')))
    return df

def to_shp(
    out_version : str, # eg v0
    SWORD_version : str,
    root_dir : Union[str,pathlib.Path] = ...,
    out_dir : Union[str,pathlib.Path] = ".",
    max_workers : int = 24,
    overwrite = False,
    ):

    if not isinstance(out_dir, pathlib.Path):
        out_dir = pathlib.Path(out_dir)
    assert out_dir.exists(), f'Out dir "{str(root_dir)}" not found!'

    shp_file = (out_dir / 'shp' / ('IRIS_shp_' + out_version + '.shp'))

    if overwrite and shp_file.exists():
        shp_file.unlink(missing_ok=True)
    elif not overwrite and shp_file.exists():
        print(f'Shapefile {str(shp_file)} exsts and overwrite flag not set. Doing nothing...')
        return None

    pandas_file = (out_dir / ('IRIS_pandas_' + out_version + '.pkl'))
    if pandas_file.exists():
        df = pickle.load(open(pandas_file,'rb'))
    else:
        df = to_df(
            out_version=out_version,
            SWORD_version=SWORD_version,
            root_dir=root_dir,
            out_dir=out_dir,
            max_workers=max_workers
            )

    sword = SWORD_IF()
    root = sword.sword_root / 'netcdf'
    dfs = []
    df.loc[:,'processed'] = True
    for file in root.iterdir():
        if '.nc' in file.name:
            dset = sword.get_dset(file)
            _df = pd.DataFrame({
                'reach_id' : dset['reaches']['reach_id'][:],
                'facc' : dset['reaches']['facc'][:],
                'width' : dset['reaches']['width'][:],
                'swordslope' : dset['reaches']['slope'][:]*1000.,
                'xmax' : dset['reaches']['x_max'][:],
                'ymax' : dset['reaches']['y_max'][:],
                'xmin' : dset['reaches']['x_min'][:],
                'ymin' : dset['reaches']['y_min'][:],
                'dist_out' : dset['reaches']['dist_out'][:],
                'wse' : dset['reaches']['wse'][:],
            })
            dset.close()
            dfs.append(_df)
    all_df = pd.concat(dfs)
    all_df.loc[:,'lon'] = all_df[['xmin', 'xmax']].mean(axis=1)
    all_df.loc[:,'lat'] = all_df[['ymin', 'ymax']].mean(axis=1)
    all_df.drop(['xmax','xmin','ymax','ymin'],axis=1,inplace=True)

    df = pd.merge(all_df,df.drop(labels=['lon','lat'],axis=1),on='reach_id',how='outer')
    df.fillna({'processed':False},inplace=True)
    geom = df.apply(lambda x: Point(x.lon,x.lat),1)
    gdf = gpd.GeoDataFrame(df,geometry=geom,crs=4326)
    
    gdf.rename(
        columns={                   #0123456789
            "across_flag":          "acr_flag",
            "along_flag":           "alg_flag",
            "combined_flag":        "com_flag",
            "avg_across_slope":     "avg_acr",
            "avg_along_slope":      "avg_alg",
            "avg_combined_slope":   "avg_com",
            "min_across_slope":     "min_acr",
            "min_along_slope":      "min_alg",
            "min_combined_slope":   "min_com",
            "max_across_slope":     "max_acr",
            "max_along_slope":      "max_alg",
            "max_combined_slope":   "max_com",
            "std_across_slope":     "std_acr",
            "std_along_slope":      "std_alg",
            "std_combined_slope":   "std_com",
            "n_across_slope":       "n_acr",
            "n_along_slope":        "n_alg",
            "n_combined_slope":     "n_com",
            "min_date_across_slope":"min_acr_dt",
            "min_date_along_slope": "min_alg_dt",
            "min_date_combined_slope":"min_com_dt",
            "max_date_across_slope":"max_acr_dt",
            "max_date_along_slope": "max_alg_dt",
            "max_date_combined_slope":"max_com_dt",
            },inplace=True)
    
    for field in ['min_acr_dt','min_alg_dt','min_com_dt','max_acr_dt','max_alg_dt','max_com_dt']:
        gdf.loc[:,field] = df.loc[:,field].apply(lambda x: x.strftime('%Y-%m-%d') if isinstance(x,date) else x)

    gdf.to_file(str(shp_file))
    return gdf
    

def to_nc(
    out_version : str, # eg v0
    SWORD_version : str,
    root_dir : Union[str,pathlib.Path] = ...,
    out_dir : Union[str,pathlib.Path] = ".",
    max_workers : int = 24,
    overwrite = False,
    ):

    assert SWORD_version in SWORD_DOI_LOOKUP.keys(), f'SWORD version {SWORD_version} not in SWORD_DOI_LOOKUP. Create new entry or use one of: {", ".join(list(SWORD_DOI_LOOKUP.keys()))}'

    if not isinstance(out_dir, pathlib.Path):
        out_dir = pathlib.Path(out_dir)
    assert out_dir.exists(), f'Out dir "{str(root_dir)}" not found!'

    nc_file = (out_dir / ('IRIS_netcdf_' + out_version + '.nc'))
    if overwrite and nc_file.exists():
        nc_file.unlink(missing_ok=True)
    elif not overwrite and nc_file.exists():
        print(f'NetCDF file {str(nc_file)} exsts and overrite flag not set. Doing nothing...')
        return None

    pandas_file = (out_dir / ('IRIS_pandas_' + out_version + '.pkl'))
    if pandas_file.exists():
        df = pickle.load(open(pandas_file,'rb'))
    else:
        df = to_df(
            out_version=out_version,
            SWORD_version=SWORD_version,
            root_dir=root_dir,
            out_dir=out_dir,
            max_workers=max_workers
            )
    
    ds = netCDF4.Dataset(nc_file, 'w',format="NETCDF4")
    ds.title = 'IRIS: ICESat-2 River Surface Slope'
    ds.subtitle = f'Supplementary data for the SWOT River Database (SWORD) {SWORD_version} ({SWORD_DOI_LOOKUP[SWORD_version]})'
    ds.history = 'Created ' + time.ctime(time.time())
    ds.institution = "DGFI-TUM"
    ds.source = 'ICESat-2 spaceborne lidar'
    ds.author = 'Daniel Scherer, Christian Schwatke, Denise Dettmering, and Florian Seitz'
    ds.contact = 'daniel.scherer@tum.de'
    ds.reference = 'Scherer D., Schwatke C., Dettmering D., and Seitz F.: ICESat-2 based River Surface Slope and its Impact on Water Level Time Series from Satellite Altimetry (doi: 10.1029/2022WR032842)'
    ds.version = out_version
    print(ds)

    c, dtype = 'reach_id', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'SWORD reach identifier'
    var[:] = df.reach_id

    c, dtype = 'lon', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Approx. centroid longitude of SWORD reach'
    var.unit = 'degrees'
    var[:] = df.lon

    c, dtype = 'lat', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Approx. centroid latitude of SWORD reach'
    var.unit = 'degrees'
    var[:] = df.lat

    c, dtype = 'across_flag', 'i1'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Flag indicating whether ICESat-2 across slope is available (1) for the reach or not (0)'
    var[:] = var[:] = df.loc[:,c].astype(int)

    c, dtype = 'along_flag', 'i1'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Flag indicating whether ICESat-2 along slope is available (1) for the reach or not (0)'
    var[:] = var[:] = df.loc[:,c].astype(int)

    c, dtype = 'combined_flag', 'i1'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Flag indicating whether ICESat-2 combined slope is available (1) for the reach or not (0)'
    var[:] = var[:] = df.loc[:,c].astype(int)

    c, dtype = 'avg_across_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Average ICESat-2 across slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'avg_along_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Average ICESat-2 along slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'avg_combined_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Average ICESat-2 combined slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'min_across_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Minimum ICESat-2 across slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'min_along_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Minimum ICESat-2 along slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'min_combined_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Minimum ICESat-2 combined slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'max_across_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Maximum ICESat-2 across slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'max_along_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Maximum ICESat-2 along slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'max_combined_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Maximum ICESat-2 combined slope for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'std_across_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'ICESat-2 across slope standard deviation for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'std_along_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'ICESat-2 along slope standard deviation for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'std_combined_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'ICESat-2 combined slope standard deviation for the reach'
    var.unit = 'mm/km'
    var[:] = df.loc[:,c].round()

    c, dtype = 'n_across_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Number of days with ICESat-2 across slope observations for the reach'
    var[:] = df.loc[:,c]

    c, dtype = 'n_along_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Number of days with ICESat-2 along slope observations for the reach'
    var[:] = df.loc[:,c].astype(int)

    c, dtype = 'n_combined_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Number of days with ICESat-2 combined slope observations for the reach'
    var[:] = df.loc[:,c]

    c, dtype = 'min_date_across_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'First date of ICESat-2 across slope observations for the reach'
    var[:] = df.loc[:,c].apply(lambda x: (x - date(2000,1,1)).days if isinstance(x,date) else np.nan)
    var.unit = 'days since 2000-01-01'

    c, dtype = 'min_date_along_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'First date of ICESat-2 along slope observations for the reach'
    var[:] = df.loc[:,c].apply(lambda x: (x - date(2000,1,1)).days if isinstance(x,date) else np.nan)
    var.unit = 'days since 2000-01-01'

    c, dtype = 'min_date_combined_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'First date of ICESat-2 combined slope observations for the reach'
    var[:] = df.loc[:,c].apply(lambda x: (x - date(2000,1,1)).days if isinstance(x,date) else np.nan)
    var.unit = 'days since 2000-01-01'

    c, dtype = 'max_date_across_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Latest date of ICESat-2 across slope observations for the reach'
    var[:] = df.loc[:,c].apply(lambda x: (x - date(2000,1,1)).days if isinstance(x,date) else np.nan)
    var.unit = 'days since 2000-01-01'

    c, dtype = 'max_date_along_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Latest date of ICESat-2 along slope observations for the reach'
    var[:] = df.loc[:,c].apply(lambda x: (x - date(2000,1,1)).days if isinstance(x,date) else np.nan)
    var.unit = 'days since 2000-01-01'

    c, dtype = 'max_date_combined_slope', 'f8'
    dim = ds.createDimension(c, None)
    var = ds.createVariable(c, dtype, (c))
    var.description = 'Latest date of ICESat-2 combined slope observations for the reach'
    var[:] = df.loc[:,c].apply(lambda x: (x - date(2000,1,1)).days if isinstance(x,date) else np.nan)
    var.unit = 'days since 2000-01-01'

    ds.close()