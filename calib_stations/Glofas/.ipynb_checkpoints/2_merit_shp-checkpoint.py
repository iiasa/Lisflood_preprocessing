"""
# -------------------------------------------------------------------------
# Name:        Use MERIT coordinates of upstream area to create shapefiles
# Purpose:     uses upstream area of MERIT (UPA) and GRDC station data
#              to create shapefiles  from rivernetwork MERIT data
#
# Author:      Peter Burek, Jesús Casado Rodríguez
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022
#
# Updated:     12/08/2024
#
# input:  grdc_MERIT_1.txt  station with new location fitted to merit UPA
# output: grdc_MERIT_2.txt: station with new location fitted to merit UPA and shapefile
#       shapefiles of basins in : ashapex_merit  e.g. grdc_basin_merit_1104800.shp

No: Number from 1 ...
GRDC_No: GRDC number
lat: original latitude from GRDC metafile
lon: original longituted from GRDC metafile
area; provided basin area from GRDC metafile
newlat: corrected latitude based on MERIT UPA dataset
newlon: corrected longitute based on MERIT UPA dataset
newarea: basin area based on MERIT UPA dataset
shapearea: area of the shape calculated with geopandas directly from shape
Lehner_shape_avail: Is there a shapefile from Lehner, 2012 to compare with?
lehner_shapearea: area of the shape (lehner, 2012) calculated with geopandas directly from shape
Indicator: indicator of similarity calculated by intersection / union of newshae and lehner shape - closer to 1 is more similar

"""

import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import pandas as pd
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


# CONFIGURATION

# input
STATION_FILE = Path('stations.csv')
LDD_FINE_FILE = Path('../data/danube_fd.tif')

# output folder
SHAPE_FOLDER = Path('./shapefiles/')


# READ INPUT DATA

# read local drainage direction map
ldd_fine = rioxarray.open_rasterio(LDD_FINE_FILE).squeeze(dim='band')
logger.info(f'Map of local drainage directions corretly read: {LDD_FINE_FILE}')

# create river network
fdir_fine = pyflwdir.from_array(ldd_fine.data,
                                ftype='d8', 
                                transform=ldd_fine.rio.transform(),
                                check_ftype=False,
                                latlon=True)

# resolution of the input map
cellsize = np.mean(np.diff(ldd_fine.x)) # degrees
cellsize_arcsec = int(np.round(cellsize * 3600, 0)) # arcsec
suffix_fine = f'{cellsize_arcsec}sec'
logger.info(f'Fine resolution is {cellsize_arcsec} arcseconds')

# read stations text file
stations = pd.read_csv(f'{STATION_FILE.stem}_{suffix_fine}.csv', index_col='ID')
logger.info(f'Table of stations correctly read: {STATION_FILE}')

# output path
SHAPE_FOLDER_FINE = SHAPE_FOLDER / suffix_fine
SHAPE_FOLDER_FINE.mkdir(parents=True, exist_ok=True)


# PROCESSING

for ID, attrs in tqdm(stations.iterrows(), total=stations.shape[0], desc='stations'):
       
    # corrected coordinates
    lat, lon = attrs[[f'lat_{suffix_fine}', f'lon_{suffix_fine}']]
    
    # boolean map of the catchment associated to the corrected coordinates
    logger.debug("basin")
    try:
        basin_arr = fdir_fine.basins(xy=(lon, lat)).astype(np.int32)
    except Exception as e:
        logger.error(f'Conversion to basin not working in catchment {ID}: {e}')
        continue

    # vectorize the boolean map into geopandas
    logger.debug("vectorize")
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
    logger.info(f'Catchment {ID} exported as shapefile: {output_file}')

