# Lisflood preprocessing

This repository contains tools useful for pre-processing the inputs of the hydrological model [LISFLOOD-OS](https://ec-jrc.github.io/lisflood/). It is part of the outcomes of the European project SEED-FD (Project - HORIZON-CL4-2023-SPACE-01).

## Installation

Get a local copy of the repository. You can either download it from GitHub or clone it with Git:

```git clone https://github.com/casadoj/lisflood-preprocessing.git```

Move to the root directory of the repository you've just copied:

```cd <YOUR_PATH>/lisflood-preprocessing/```

Install the package with PiP:

```pip install .```

## Tools

### [`lfcoords`](./src/lisfloodpreprocessing/lfcoords.py)

This tool finds the appropriate coordinates in the LISFLOOD river network of any point, provided that the catchment area is known. A thourough explanation of the method can be found in [Burek and Smilovic (2023)](https://essd.copernicus.org/articles/15/5617/2023/).

First, it uses the original coordinates and catchment area to find the most accurate pixel in a high-resolution map. [Burek and Smilovic (2023)](https://essd.copernicus.org/articles/15/5617/2023/) use MERIT [(Yamazaki et al., 2019)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019WR024873), which has a spatial resolution of 3 arc-seconds. The result of this first step is, for every point, a new value of coordinates and area, and a shapefile of the catchment polygon in high-resolution.

Second, it finds the pixel in the low-resolution grid (LISFLOOD static maps) that better matches the catchment shape derived in the previous step. As a result, for each point we obtained a new value of coordinates and area, and a new shapefile of the catchment polygon in low-resolution.

Reservoirs in LISFLOOD need to be located one pixel downstream of their actual location. This is because LISFLOOD reads the reservoir inflow from one pixel upstream of the location defined in the reservoir static map, and releases the outflow at the reservoir pixel. The tool includes a flag (`-r`) to specify that the input points are reservoirs; if that is flagged, the coordinates in the CSV output file are corrected to refer to one pixel downstream, i.e., they can be directly used in LISFLOOD.

#### Usage

Once the package is installed, the tool can be executed from the command prompt by indicating a configuration file:

```bash
lfcoords --config-file config.yml
```

If the input points are reservoirs, you can use the flag `-r` to get the coordinates one pixel downstream of the original result:

```bash
lfcoords --config-file config.yml -r
```

##### Configuration file

The configuration file defines the input files, the folder where the resulting shapefiles will be saved, and some thresholds used in the process. A template of the configuration file can be found [here](./src/lisfloodpreprocessing/config.yml). Below you find an example:

```yml
input:
    points: stations.csv # ID, lat, lon, area
    ldd_fine: ../data/danube_fd.tif
    upstream_fine: ../data/ups_danube_3sec.tif # km2
    ldd_coarse: ldd_3min.tif
    upstream_coarse: upArea_repaired.nc # m2
            
output_folder: ./shapefiles/

conditions:
    min_area: 100 # km2
    abs_error: 50 # km2
    pct_error: 1 # %
```

##### Inputs

The tool requires 5 inputs:

* A CSV file of the stations to be located in the LISFLOOD grid. This file must contain four columns with four specific names: 'ID', 'area' in km2, 'lat', 'lon'. Below you find an example of the stations CSV file:

```csv
ID,area,lat,lon
429,35399,49.018,12.144
436,26448,48.947,12.015
439,37687,48.88,12.747
```

* A map of the local drainage directions in high-resolution, e.g., MERIT.
* A map of the upstream area in high-resolution, e.g., MERIT. The units of this map must be km2, same units as the _area_ field in the CSV file.
* A map of the local drainage directions in low-resolution, i.e., the LISFLOOD static map.
* A map of the upstream area in low-resolution, i.e., the LISFLOOD static map. The units of this map are m2 (instead of km2), as these are the units used in LISFLOOD; the code converts internally this map into km2.

All maps can be provided either in TIFF or NetCDF format.

##### Outputs

The main output is a new **CSV file** saved in the same directory as the input CSV file and named similarly, but with a suffix indicating the resolution of the LISFLOOD grid. For instance, in the configuration file above the input CSV file is named _stations.csv_ and the resolution of the LISFLOOD grid is 3 arcmin, so the output CSV file will be named _stations_3min.csv_. The CSV contains 6 new columns defining the coordinates and catchment area in both the high-resolution (`3sec` in the example) and low-resolution grids (`3min` in the example). Example:

```csv
ID,area,area_3min,area_3sec,lat,lat_3min,lat_3sec,lon,lon_3min,lon_3sec
429,35399,35216,35344,49.018,49.025,49.022083,12.144,12.125,12.14375
436,26448,26334,26394,48.947,48.925,48.94625,12.015,12.025,12.014583
439,37687,37540,37605,48.88,48.925,48.879583,12.747,12.675,12.74625
```

Besides, the tool creates **shapefiles** of the catchment polygons derived for both the high and low resolution grids. The shapefiles are saved in two subdirectories inside the `output_folder` directory defined in the configuration file. In each of these subdirectories, there will be one file for each station.