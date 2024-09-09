import numpy as np
import pandas as pd
import rioxarray
import pyflwdir
from tqdm.auto import tqdm
from pathlib import Path

import warnings
warnings.filterwarnings("ignore")

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
logger = logging.getLogger(__name__)

from utils import find_pixel, catchment_polygon 
from __init__ import Config


# CONFIGURATION

cfg = Config('config.yml')


# READ INPUT DATA

# read upstream map with fine resolution
upstream_fine = rioxarray.open_rasterio(cfg.UPSTREAM_FINE).squeeze(dim='band')
logger.info(f'Map of upstream area corretly read: {cfg.UPSTREAM_FINE}')

# read local drainage direction map
ldd_fine = rioxarray.open_rasterio(cfg.LDD_FINE).squeeze(dim='band')
logger.info(f'Map of local drainage directions corretly read: {cfg.LDD_FINE}')

# resolution of the input map
cellsize = np.mean(np.diff(upstream_fine.x)) # degrees
cellsize_arcsec = int(np.round(cellsize * 3600, 0)) # arcsec
suffix_fine = f'{cellsize_arcsec}sec'
logger.info(f'The resolution of the finer grid is {cellsize_arcsec} arcseconds')

# read stations text file
stations = pd.read_csv(cfg.STATIONS, index_col='ID')
logger.info(f'Table of stations correctly read: {cfg.STATIONS}')


# PROCESSING

# add columns to the table of stations
new_cols = [f'{col}_{suffix_fine}' for col in stations.columns]
stations[new_cols] = np.nan

# create river network
fdir_fine = pyflwdir.from_array(ldd_fine.data,
                                ftype='d8', 
                                transform=ldd_fine.rio.transform(),
                                check_ftype=False,
                                latlon=True)

# output path
SHAPE_FOLDER_FINE = cfg.SHAPE_FOLDER / suffix_fine
SHAPE_FOLDER_FINE.mkdir(parents=True, exist_ok=True)

for ID, attrs in tqdm(stations.iterrows(), total=stations.shape[0], desc='stations'):  

    # reference coordinates and upstream area
    lat_ref, lon_ref, area_ref = attrs[['lat', 'lon', 'area']]

    # search new coordinates in an increasing range
    ranges = [55, 101, 151]
    penalties = [500, 500, 1000]
    factors = [2, .5, .25]
    acceptable_errors = [50, 80, np.nan]
    for rangexy, penalty, factor, max_error in zip(ranges, penalties, factors, acceptable_errors):
        logger.debug(f'Set range to {rangexy}')
        lat, lon, error = find_pixel(upstream_fine, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=penalty, factor=factor)
        if error <= max_error:
            break

    # update new columns in 'stations'
    stations.loc[ID, new_cols] = [round(lat, 6), round(lon, 6), int(upstream_fine.sel(y=lat, x=lon).item())]
    
    # boolean map of the catchment associated to the corrected coordinates
    basin_arr = fdir_fine.basins(xy=(lon, lat)).astype(np.int32)

    # vectorize the boolean map into geopandas
    basin_gdf = catchment_polygon(basin_arr.astype(np.int32),
                                  transform=ldd_fine.rio.transform(),
                                  crs=ldd_fine.rio.crs,
                                  name='ID')
    basin_gdf['ID'] = ID
    basin_gdf.set_index('ID', inplace=True)
    basin_gdf[attrs.index] = attrs.values

    # export shape file
    output_file = SHAPE_FOLDER_FINE / f'{ID}.shp'
    basin_gdf.to_file(output_file)
    logger.info(f'\nCatchment {ID} exported as shapefile: {output_file}')

# export results
stations.sort_index(axis=1, inplace=True)
output_csv = f'{cfg.STATIONS.stem}_{suffix_fine}.csv'
stations.to_csv(output_csv)
logger.info(f'Coordinates an upstream area in the finer grid have been exported to: {output_csv}')