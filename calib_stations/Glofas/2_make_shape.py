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


# CONFIGURATION

# input
STATION_FILE = 'stations_MERIT_xarray.csv' # output of step 1
LDD_FILE = '../data/danube_fd.tif'

# output
SHAPE_FOLDER = Path('../shape_glofas_3sec/')


# READ INPUT DATA

# read stations text file
stations = pd.read_csv(STATION_FILE, index_col='ID')
logger.info(f'Table of stations correctly read: {STATION_FILE}')

# read local drainage direction map
ldd = rioxarray.open_rasterio(LDD_FILE).squeeze(dim='band')
ldd = ldd.rename({'x': 'lon', 'y': 'lat'})
logger.info(f'Map of local drainage directions corretly read: {LDD_FILE}')

# parse the LDD map to be usable later on
logger.debug("make ldd")
flw = pyflwdir.from_array(ldd.data, #flwdir, 
                          ftype='d8', 
                          transform=ldd.rio.transform(),
                          check_ftype=False,
                          latlon=True)
logger.debug("done ldd")


# PROCESSING

for ID, attrs in tqdm(stations.iterrows(), total=stations.shape[0], desc='stations'):
       
    # corrected coordinates
    lat, lon = attrs[['lat_new', 'lon_new']]
    
    # boolean map of the catchment associated to the corrected coordinates
    logger.debug("basin")
    try:
        basin_arr = flw.basins(xy=(lon, lat)).astype(np.int32)
    except Exception as e:
        logger.error(f'Conversion to basin not working in catchment {ID}: {e}')
        continue

    # vectorize the boolean map into geopandas
    logger.debug("vectorize")
    basin_gdf = catchment_polygon(basin_arr.astype(np.int32),
                                  transform=ldd.rio.transform(),
                                  crs=ldd.rio.crs,
                                  name='ID')
    basin_gdf['ID'] = ID
    basin_gdf.set_index('ID', inplace=True)
    basin_gdf[attrs.index] = attrs.values

    # export shape file
    output_file = SHAPE_FOLDER / f'{ID}.shp'
    basin_gdf.to_file(output_file)
    logger.info(f'Catchment {ID} exported as shapefile: {output_file}')