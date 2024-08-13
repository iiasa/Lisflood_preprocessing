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


# CONFIGURATION

# input
STATION_FILE = Path('stations.csv')
UPSTREAM_FINE_FILE = Path('../data/ups_danube_3sec.tif') # MERIT Yamazaki et al 2019 
LDD_FINE_FILE = Path('../data/danube_fd.tif')

# output
SHAPE_FOLDER = Path('./shapefiles/')


# READ INPUT DATA

# read upstream map with fine resolution
upstream_fine = rioxarray.open_rasterio(UPSTREAM_FINE_FILE).squeeze(dim='band')
logger.info(f'Map of upstream area corretly read: {UPSTREAM_FINE_FILE}')

# read local drainage direction map
ldd_fine = rioxarray.open_rasterio(LDD_FINE_FILE).squeeze(dim='band')
logger.info(f'Map of local drainage directions corretly read: {LDD_FINE_FILE}')

# resolution of the input map
cellsize = np.mean(np.diff(upstream_fine.x)) # degrees
cellsize_arcsec = int(np.round(cellsize * 3600, 0)) # arcsec
suffix_fine = f'{cellsize_arcsec}sec'
logger.info(f'The resolution of the finer grid is {cellsize_arcsec} arcseconds')

# read stations text file
stations = pd.read_csv(STATION_FILE, index_col='ID')
logger.info(f'Table of stations correctly read: {STATION_FILE}')


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
SHAPE_FOLDER_FINE = SHAPE_FOLDER / suffix_fine
SHAPE_FOLDER_FINE.mkdir(parents=True, exist_ok=True)

for ID, attrs in tqdm(stations.iterrows(), total=stations.shape[0], desc='stations'):  

    # reference coordinates and upstream area
    lat_ref, lon_ref, area_ref = attrs[['lat', 'lon', 'area']]

    # search range in cells: 55 = around 5km
    rangexy = 55
    logger.debug(f'Set range to {rangexy}')
    lat, lon, error = find_pixel(upstream_fine, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=500, factor=2)
    
    # if still big error, increase range
    if error > 50:
        rangexy = 101
        logger.debug(f'Increase range to {rangexy}')
        lat, lon, error = find_pixel(upstream_fine, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=500, factor=0.5)

        # if still big error increase range
        if error > 80:
            rangexy = 151
            logger.debug(f'Increase range to {rangexy}')
            lat, lon, error = find_pixel(upstream_fine, lat_ref, lon_ref, area_ref, rangexy=rangexy, penalty=1000, factor=0.25)

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
output_csv = f'{STATION_FILE.stem}_{suffix_fine}.csv'
stations.to_csv(output_csv)
logger.info(f'Coordinates an upstream area in the finer grid have been exported to: {output_csv}')