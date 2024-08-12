# ***
# Name:        Find MERIT coordinates
# Purpose:     uses upstream area of MERIT (UPA) and GRDC station data
#              to check and correct station location
# 
"""
# Author:      Peter Burek, Jesús Casado Rodríguez
# 
# Created:     15/05/2022
# Copyright:   (c) PB 2022
#
# Updated:     12/08/2024
# 
# input:  grdc_2022_10577.txt   10577 station datasets >= 10km2 upstream area or no area provided
# output: grdc_MERIT_1.txt: station with new location fitted to merit UPA
# 
# No: Number from 1 ...
# GRDC_No: GRDC number
# lat: original latitude from GRDC metafile
# lon: original longituted from GRDC metafile
# newlat: corrected latitude based on MERIT UPA dataset
# newlon: corrected longitute based on MERIT UPA dataset
# area; provided basin area from GRDC metafile
# newarea: basin area based on MERIT UPA dataset
# UPS_Indicator:  min error in % from MERIT UPA to provided basin area
# dist_Indicator: distance to original pour point in [unit:100m]
# Indicator:  ranking criteria: UPS_Indicator + 2 x dist_indicator
# 
# ***
"""

import numpy as np
import pandas as pd
import xarray as xr
#import rasterio
import rioxarray
from tqdm.auto import tqdm
from typing import Tuple

import warnings
warnings.filterwarnings('ignore')

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
logger = logging.getLogger(__name__)

from utils import find_pixel


# CONFIGURATION

# input
STATION_FILE = 'stations.csv' # table of original coordinates (lat and lon) and catchment area in km2
UPSTREAM_FILE = "../data/ups_danube_3sec.tif" # MERIT Yamazaki et al 2019 - upstream area in km2

# output
OUTPUT_FILE = 'stations_MERIT_xarray.csv'


# READ INPUT DATA

# read stations text file
stations = pd.read_csv(STATION_FILE, index_col='ID')
new_cols = [f'{col}_new' for col in stations.columns]
stations[new_cols] = np.nan
logger.info(f'Table of stations correctly read: {STATION_FILE}')

# read upstream map
upstream = rioxarray.open_rasterio(UPSTREAM_FILE).squeeze(dim='band')
upstream = upstream.rename({'x': 'lon', 'y': 'lat'})
logger.info(f'Map of upstream area corretly read: {UPSTREAM_FILE}')


# PROCESSING

for ID, attrs in tqdm(stations.iterrows(), total=stations.shape[0], desc='stations'):  

    # reference coordinates and upstream area
    lat_ref, lon_ref, area_ref = attrs[['lat', 'lon', 'area']]

    # search range in cells: 55 = around 5km
    rangexy = 55
    logger.debug(f'Set range to {rangexy}')
    lat, lon, error = find_pixel(upstream, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=500, factor=2)
    
    # if still big error, increase range
    if error > 50:
        rangexy = 101
        logger.debug(f'Increase range to {rangexy}')
        lat, lon, error = find_pixel(upstream, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=500, factor=0.5)

    # if still big error increase range
        if error > 80:
            rangexy = 151
            logger.debug(f'Increase range to {rangexy}')
            lat, lon, error = find_pixel(upstream, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=1000, factor=0.25)

    # find new coordinates and its associated upstream area
    stations.loc[ID, ['lat_new', 'lon_new', 'area_new']] = [round(lat, 6), round(lon, 6), int(upstream.sel(lat=lat, lon=lon).item())]

# export results
stations.sort_index(axis=1, inplace=True)
stations.to_csv(OUTPUT_FILE)
logger.info(f'Results have been exported to: {OUTPUT_FILE}')