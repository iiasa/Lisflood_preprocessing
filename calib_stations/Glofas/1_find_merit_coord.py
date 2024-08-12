#!/usr/bin/env python
# coding: utf-8

# ***
# Name:        Find MERIT coordinates
# Purpose:     uses upstream area of MERIT (UPA) and GRDC station data
#              to check and correct station location
# 
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


def find_pixel(
    upstream: xr.DataArray,
    lat: int,
    lon: int,
    area: float,
    rangexy: int = 55,
    penalty: int = 500,
    factor: int = 2,
    distance_scaler: float = .92,
    error_threshold: int = 50
) -> Tuple:
    """
    Find the coordinates of the pixel in the upstream map with a smaller error compared with a reference area.
    
    Parameters:
    -----------
    upstream: xr.DataArray
        The upstream data containing latitude and longitude coordinates.
    lat: float
        The original latitude value.
    lon: float
        The original longitude value.
    area: float
        The reference area to calculate percent error.
    rangexy: int, optional
        The range in both x and y directions to search for the new location.
    penalty: int, optional
        The penalty value to add to the distance when the percent error is too high.
    error_threshold: float, optional
        The threshold for the percent error to apply the penalty.
    factor: int, optional
        The factor to multiply with the distance for the error calculation.
    distance_scaler: float, optional
        The scaling factor for the distance calculation in pixels.
    
    Returns:
    --------
    lat_new : float
        The latitude of the new location.
    lon_new : float
        The longitude of the new location.
    min_error : float
        The minimum error value at the new location.
    """

    # find coordinates of the nearest pixel in the map
    nearest_pixel = upstream.sel(lat=lat, lon=lon, method='nearest')
    lat_orig, lon_orig = [nearest_pixel[coord].item() for coord in ['lat', 'lon']]

    # extract subset of the upstream map
    cellsize = np.mean(np.diff(upstream.lon.data))
    delta = rangexy * cellsize + 1e-6
    upstream_sel = upstream.sel(lat=slice(lat_orig + delta, lat_orig - delta),
                                lon=slice(lon_orig - delta, lon_orig + delta))

    # percent error in catchment area
    error = 100 * (1 - upstream_sel / area)

    # distance from the original pixel (in pixels)
    i = np.arange(-rangexy, rangexy + 1)
    ii, jj = np.meshgrid(i, i)
    distance = xr.DataArray(data=np.sqrt(ii**2 + jj**2) * distance_scaler, coords=upstream_sel.coords, dims=upstream_sel.dims)
    # penalise if error is too big
    distance = distance.where(error <= error_threshold, distance + penalty)

    # update error based on distance
    error += factor * distance

    # coordinates of the new location
    min_error = error.where(error == error.min(), drop=True)
    lat_new, lon_new = [min_error[coord].item() for coord in ['lat', 'lon']]
    
    return lat_new, lon_new, min_error.item()


# CONFIGURACIÓN

# input

STATION_FILE = 'stations.csv' # table of original coordinates (lat and lon) and catchment area in km2
MAP = "../data/ups_danube_3sec.tif" # MERIT Yamazaki et al 2019 - upstream area in km2

# output
OUTPUT_FILE = 'stations_Merit_xarray.csv'


# READ INPUT DATA

# read stations text file
stations = pd.read_csv(STATION_FILE, index_col='ID')
new_cols = [f'{col}_new' for col in stations.columns]
stations[new_cols] = np.nan
logger.info(f'Table of stations correctly read: {STATION_FILE}')

# read upstream map
upstream = rioxarray.open_rasterio(MAP).squeeze(dim='band')
upstream = upstream.rename({'x': 'lon', 'y': 'lat'})
logger.info(f'Map or upstream area corretly read: {STATION_FILE}')


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