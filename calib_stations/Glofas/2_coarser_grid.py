"""
# -------------------------------------------------------------------------
# Name:        Shapefile in low-res
# Purpose:     creates shapefiles and a list in lower resolution (5 or 30 arcmin)
#
# Author:      Peter Burek, Jesús Casado Rodríguez
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022
#
# Update:      13/08/2024
#
# input:  results/grdc_shape_allend_1.txt  station with new location fitted to merit UPA from 3_makeshape.py
# output: basins_30min.txt or basins_5min.txt: station with new location fitted to 30 arcmin or 5 arcmin
#          shapefiles of low-res basins in : ashape30min or ashape5min e.g. grdc_basin_30min_basin_1104150.shp

basins_30min.txt
----------------
No: Number from 0 ...
GRDC_No: GRDC number
similarity: similarity of 30min low-res shape with high-res shape from 3arcsec
areaGRDC: area provided by GRDC
area: area from high-res UPA MERIT at pour point
lat: coorected latitude on high-res
lon: corrected longitude on high-res
area30min : area on pour point of best fitting 30arcmin grid cell
lat30min: latitude on 30arcmin
lon30min: longitude on 30arcmin
lat30move: latitude moved to grid cell centre
lon30move: longitude moved to grid cell centre
indloc: locatation of indice which is taken as nbest fit: a 5x5 gridcell surrounding is taken numbered 0-24
upsloc: location (index of 0-24) of best fitting upsstream area
ind:    scoring result for best upstream fitting based on upstream area fitting and shape similarity
ratio: scoring on upstream area fitting between low-res area and high-res area
ups1:   low-res upstream area of best upstream fitting
indshape: scoring on similarity between low-res shape and high res shape
shapeloc: location (index of 0-24) of best fitting shape similarity
ind:    scoring result for best upstream fitting based on upstream area fitting and shape similarity
ratio: scoring on upstream area fitting between low-res area and high-res area
ups1:   low-res upstream area of best shape fitting
indshape: scoring on similarity between low-res shape and high res shape

# uses
# ldd_30min.tif: river network in pcraster format
30min from Döll, P., Lehner, B. (2002): Validation of a new global 30-min drainage direction map. Journal of Hydrology, 258(1-4), 214-23
5min ldd from: Eilander, D., van Verseveld, W., Yamazaki, D., Weerts, A., Winsemius, H. C., and Ward, P. J.: A hydrography upscaling method for scale-invariant parametrization of distributed hydrological models, Hydrol. Earth Syst. Sci., 25, 5287-5313, 10.5194/hess-25-5287-2021, 2021.
# ups_30min.tif: upstream area in [km]

# Shapefiles (smooth) with high resolution from: /ashape2_merit/grdc_basin_merit_


# ----------------------------------------------------------------------
"""

import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray
import pyflwdir
from pathlib import Path
from tqdm.auto import tqdm

import warnings
warnings.filterwarnings("ignore")

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
logger = logging.getLogger(__name__)

from utils import catchment_polygon
from __init__ import Config


### CONFIGURATION

cfg = Config('config.yml')

# input data
STATION_FILE = 'stations_MERIT.csv' # output of step 1


### READ INPUT DATA

# read upstream area map of coarse grid
upstream_coarse = rioxarray.open_rasterio(cfg.UPSTREAM_COARSE).squeeze(dim='band')
logger.info(f'Map of upstream area corretly read: {cfg.UPSTREAM_COARSE}')

# read local drainage direction map
ldd_coarse = rioxarray.open_rasterio(cfg.LDD_COARSE).squeeze(dim='band')
logger.info(f'Map of local drainage directions correctly read: {cfg.LDD_COARSE}')

# create river network
fdir_coarse = pyflwdir.from_array(ldd_coarse.data,
                          ftype='ldd',
                          transform=ldd_coarse.rio.transform(),
                          check_ftype=False,
                          latlon=True)

# boundaries of the input maps
lon_min, lat_min, lon_max, lat_max = np.round(ldd_coarse.rio.bounds(), 6)

# resolution of the input maps
cellsize = np.round(np.mean(np.diff(ldd_coarse.x)), 6) # degrees
cellsize_arcmin = int(np.round(cellsize * 60, 0)) # arcmin
suffix_coarse = f'{cellsize_arcmin}min'
logger.info(f'Coarse resolution is {cellsize_arcmin} arcminutes')

# read stations text file
stations = pd.read_csv(STATION_FILE, index_col='ID')
logger.info(f'Table of stations correctly read: {STATION_FILE}')

# extract resolution of the finer grid
suffixes = [col.split('_')[1] for col in stations.columns if '_' in col]
suffix_fine = sorted(set(suffixes))[0]
cols_fine = [f'{col}_{suffix_fine}' for col in ['lat', 'lon', 'area']]

# add new columns
cols_coarse = [f'{col}_{suffix_coarse}' for col in ['lat', 'lon', 'area']]
stations[cols_coarse] = np.nan

# output folders
SHAPE_FOLDER_FINE = cfg.SHAPE_FOLDER / suffix_fine
SHAPE_FOLDER_COARSE = cfg.SHAPE_FOLDER / suffix_coarse
for folder in [SHAPE_FOLDER_FINE, SHAPE_FOLDER_COARSE]:
    folder.mkdir(parents=True, exist_ok=True)


### PROCESSING

# search range of 5x5 array -> this is where the best point can be found in the coarse grid
rangexy = np.linspace(-2, 2, 5) * cellsize # arcmin

for ID, attrs in tqdm(stations.iterrows(), total=stations.shape[0], desc='stations'):
    
    # real upstream area
    area_ref = attrs['area']
    
    # coordinates and upstream area in the fine grid
    lat_fine, lon_fine, area_fine = attrs[[f'{col}_{suffix_fine}' for col in ['lat', 'lon', 'area']]]
       
    if (area_ref < cfg.MIN_AREA) or (area_fine < cfg.MIN_AREA):
        logger.info(f'The catchment area of station {ID} is smaller than the minimum of {cfg.MIN_AREA} km2')
        continue
                        
    # import shapefile of catchment polygon
    shapefile = SHAPE_FOLDER_FINE / f'{ID}.shp'
    try:
        basin_fine = gpd.read_file(shapefile)
        logger.info(f'Catchment polygon correctly read: {shapefile}')
    except OSError as e:
        logger.error(f'Error reading {shapefile}: {e}')
        continue
    except Exception as e:  # This will catch other exceptions that might occur.
        logger.error(f'An unexpected error occurred while reading {shapefile}: {e}')
        continue

    # find ratio
    logger.debug('Start search')
    inter_vs_union, area_ratio, area_lisf = [], [], []
    for Δlat in rangexy:
        for Δlon in rangexy:
            lon = lon_fine + Δlon
            lat = lat_fine + Δlat
            basin = catchment_polygon(fdir_coarse.basins(xy=(lon, lat)).astype(np.int32),
                                      transform=ldd_coarse.rio.transform(), 
                                      crs=ldd_coarse.rio.crs)

            # calculate union and intersection of shapes
            intersection = gpd.overlay(basin_fine, basin, how='intersection')
            union = gpd.overlay(basin_fine, basin, how='union')
            inter_vs_union.append(intersection.area.sum() / union.area.sum())

            # get upstream area (km2) of coarse grid (LISFLOOD)
            area = upstream_coarse.sel(x=lon, y=lat, method='nearest').item() * 1e-6
            area_lisf.append(area)

            # ratio between reference and coarse area
            if area_ref == 0 or area == 0:
                ratio = 0
            else:
                ratio = area_ref / area if area_ref < area else area / area_ref
            area_ratio.append(ratio)
    logger.debug('End search')

    # maximum of shape similarity and upstream area accordance
    i_shape = np.argmax(inter_vs_union)
    area_shape = area_lisf[i_shape]
    i_centre = int(len(rangexy)**2 / 2) # middle point
    area_centre = area_lisf[i_centre]           
    # use middle point if errors are small
    abs_error = abs(area_shape - area_centre)
    pct_error = 100 * abs(1 - area_centre / area_shape)
    if (abs_error <= cfg.ABS_ERROR) and (pct_error <= cfg.PCT_ERROR):
        i_shape = i_centre
        area_shape = area_centre

    #i_ratio = np.argmax(area_ratio)          

    # coordinates in the fine resolution
    i = i_shape // len(rangexy)
    j = i_shape % len(rangexy)
    lat = lat_fine + rangexy[i]
    lon = lon_fine + rangexy[j]

    # coordinates and upstream area on coarse resolution
    area = upstream_coarse.sel(x=lon, y=lat, method='nearest')
    area_coarse = area.item() * 1e-6
    lon_coarse = area.x.item()
    lat_coarse = area.y.item()

    # derive catchment polygon from the selected coordinates
    basin_coarse = catchment_polygon(fdir_coarse.basins(xy=(lon_coarse, lat_coarse)).astype(np.int32),
                                     transform=ldd_coarse.rio.transform(), 
                                     crs=ldd_coarse.rio.crs,
                                     name='ID')
    basin_coarse['ID'] = ID
    basin_coarse.set_index('ID', inplace=True)
    basin_coarse[cols_fine] = lat_fine, lon_fine, area_fine
    basin_coarse[cols_coarse] = lat_coarse, lon_coarse, area_coarse

    # export shapefile
    output_shp = SHAPE_FOLDER_COARSE / f'{ID}.shp'
    basin_coarse.to_file(output_shp)
    logger.info(f'Catchment {ID} exported as shapefile: {output_shp}')

    # update new columns in 'stations'
    stations.loc[ID, cols_coarse] = lat_coarse, lon_coarse, area_coarse

# export results
stations.sort_index(axis=1, inplace=True)
output_csv = f'{Path(STATION_FILE).stem}_{suffix_coarse}.csv'
stations.to_csv(output_csv)
logger.info(f'Results have been exported to: {output_csv}')


# In[ ]:




