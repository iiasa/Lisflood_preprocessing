"""
# -------------------------------------------------------------------------
# Name:        Find MERIT coordinates
# Purpose:     uses upstream area of MERIT (UPA) and GRDC station data
#              to check and correct station location
#
# Author:      PB
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022

input:  grdc_2022_10577.txt   10577 station datasets >= 10km2 upstream area or no area provided
output: grdc_MERIT_1.txt: station with new location fitted to merit UPA

No: Number from 1 ...
GRDC_No: GRDC number
lat: original latitude from GRDC metafile
lon: original longituted from GRDC metafile
newlat: corrected latitude based on MERIT UPA dataset
newlon: corrected longitute based on MERIT UPA dataset
area; provided basin area from GRDC metafile
newarea: basin area based on MERIT UPA dataset
UPS_Indicator:  min error in % from MERIT UPA to provided basin area
dist_Indicator: distance to original pour point in [unit:100m]
Indicator:  ranking criteria: UPS_Indicator + 2 x dist_indicator

# ----------------------------------------------------------------------
"""
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import shapes

from osgeo import gdal
from osgeo import osr
from osgeo import gdalconst

import sys
import os


#----------------------------------------------
# INPUT
# MERIT Yamazaki et al 2019 - upstream area in km2

glofas_stations = "metastation_45.txt"
shapefolder = "../shape_glofas_3sec/shape_3sec_no_"
upsname = "../data/ups_danube_3sec.tif"

#OUTPUT
glofas_Merit = "glofas_Merit_2.txt"

# --------------------------------------------------------------------------------
# cell size: 3 arcsec
cell = 0.000833333333333333333
# 1/cell
invcell = 1200
# search range in cells: 55 = around 5km
rangexy=55

f = open(glofas_stations, "r")
glofas = f.readlines()
f.close()
header = glofas[0].rstrip()
glofas = glofas[1:]

header += "\tnewlat\tnewlon\tnewarea\n"
f = open(glofas_Merit, "w")
f.write(header)
f.close()

# -----
# load upstream
print ("read ups")
src = rasterio.open(upsname, "r")
ups= src.read(1)
transform =src.transform
latlon = src.crs.to_epsg()
crs = src.crs
src.close()
print ("done read ups")

# -----------------------------------

for stationNo in range(len(glofas)):

    station = glofas[stationNo].split("\t")
    upsreal = float(station[6])
    # upstream area from provider

    coord = [float(station[4]), float(station[5])]
    # lat lon
    glofas_no =  station[1]


    top = round(transform[5],3)

    left = round(transform[2],3)
    col = ups.shape[1]
    row = ups.shape[0]

    col1 = int((coord[1] - left) * invcell)
    row1 = int((top -coord[0]) * invcell)
    ups1 = ups[row1,col1]

    rangexy = 55
    upsups = np.zeros((rangexy*2+1,rangexy*2+1))
    ind = np.zeros((rangexy*2+1,rangexy*2+1))
    upsind = np.zeros((rangexy*2+1,rangexy*2+1))
    diffind = np.zeros((rangexy*2+1,rangexy*2+1))

    colcol = np.arange(col1-rangexy,col1+rangexy+1)
    rowrow = np.arange(row1-rangexy,row1+rangexy+1)

    j =0
    for y in rowrow:
        i = 0
        for x in colcol:
            upsind[j, i] = 100 * np.abs(1 - ups[y, x] / upsreal)
            upsups[j,i]= ups[y,x]
            diff = np.sqrt((rangexy-i)**2+(rangexy-j)**2)*0.9
            diffind[j,i] = np.sqrt((rangexy - i) ** 2 + (rangexy - j) ** 2) * 0.92
            # if upsind> 50 diff gets a penalty
            if upsind[j, i]>50:
                diffind[j, i] = diffind[j, i] + 500
            ind[j,i] = upsind[j, i] + 2 * diffind[j,i]

            i = i +1
        j = j + 1

    minxy = np.where(ind==np.min(ind))
    y=minxy[0][0]
    x=minxy[1][0]
    j = rowrow[y]
    i = colcol[x]

    yy=coord[0]+(rangexy-y)*cell
    xx=coord[1]-(rangexy-x)*cell
    ups2 = ups[j,i]


    #------------------------------------------------------
    # if still big error increase range
    if ind[y,x] > 50:
        print ("increase range")
        rangexy = 101
        upsups = np.zeros((rangexy*2+1,rangexy*2+1))
        ind = np.zeros((rangexy*2+1,rangexy*2+1))
        upsind = np.zeros((rangexy*2+1,rangexy*2+1))
        diffind = np.zeros((rangexy*2+1,rangexy*2+1))

        colcol = np.arange(col1-rangexy,col1+rangexy+1)
        rowrow = np.arange(row1-rangexy,row1+rangexy+1)

        j =0
        for y in rowrow:
            i = 0
            for x in colcol:
                upsind[j, i] = 100 * np.abs(1 - ups[y, x] / upsreal)
                upsups[j,i]= ups[y,x]
                diff = np.sqrt((rangexy-i)**2+(rangexy-j)**2)*0.9
                diffind[j,i] = np.sqrt((rangexy - i) ** 2 + (rangexy - j) ** 2) * 0.92
                # if upsind> 50 diff gets a penalty
                if upsind[j, i] > 50:
                    diffind[j, i] = diffind[j, i] + 500
                ind[j,i] = upsind[j, i] + 0.5 * diffind[j,i]

                i = i +1
            j = j + 1

        minxy = np.where(ind==np.min(ind))
        y=minxy[0][0]
        x=minxy[1][0]
        j = rowrow[y]
        i = colcol[x]

        yy=coord[0]+(rangexy-y)*cell
        xx=coord[1]-(rangexy-x)*cell
        ups2 = ups[j,i]
    #-------------------------------------------------

    # ------------------------------------------------------
    # if still big error increase range
        if ind[y, x] > 80:
            print("increase range2")
            rangexy = 151
            upsups = np.zeros((rangexy * 2 + 1, rangexy * 2 + 1))
            ind = np.zeros((rangexy * 2 + 1, rangexy * 2 + 1))
            upsind = np.zeros((rangexy * 2 + 1, rangexy * 2 + 1))
            diffind = np.zeros((rangexy * 2 + 1, rangexy * 2 + 1))

            colcol = np.arange(col1 - rangexy, col1 + rangexy + 1)
            rowrow = np.arange(row1 - rangexy, row1 + rangexy + 1)

            j = 0
            for y in rowrow:
                i = 0
                for x in colcol:
                    upsind[j, i] = 100 * np.abs(1 - ups[y, x] / upsreal)
                    upsups[j, i] = ups[y, x]
                    diff = np.sqrt((rangexy - i) ** 2 + (rangexy - j) ** 2) * 0.9
                    diffind[j, i] = np.sqrt((rangexy - i) ** 2 + (rangexy - j) ** 2) * 0.92
                    # if upsind> 50 diff gets a penalty
                    if upsind[j, i] > 50:
                        diffind[j, i] = diffind[j, i] + 1000
                    ind[j, i] = upsind[j, i] + 0.25 * diffind[j, i]

                    i = i + 1
                j = j + 1

            minxy = np.where(ind == np.min(ind))
            y = minxy[0][0]
            x = minxy[1][0]
            j = rowrow[y]
            i = colcol[x]

            yy = coord[0] + (rangexy - y) * cell
            xx = coord[1] - (rangexy - x) * cell
            ups2 = ups[j, i]
            ii =1
    # -------------------------------------------------

    #-------------------------------------------------
    s = str(stationNo)  + "\t"+ str(glofas_no) + "\t" + station[7] + "\t" + station[8] + "\t"
    s = s + f"{yy:.4f}" + "\t" + f"{xx:.4f}" + "\t" + f"{ups2:.0f}"+ "\t"
    s = s + f"{coord[0]:.4f}"+ "\t" +f"{coord[1]:.4f}" + "\t"+f"{upsreal:.0f}"
    print (s)

    #header += "\tnewlat\tnewlon\tnewarea"
    s = glofas[stationNo].rstrip()
    s = s + "\t" + f"{yy:.4f}" + "\t" + f"{xx:.4f}" + "\t" + f"{ups2:.0f}"+ "\n"
    #s = s + "\t" + str(upsind[y,x]) +"\t"+str(diffind[y,x])+"\t"+str(ind[y,x]) + "\n"
    f = open(glofas_Merit, "a")
    f.write(s)
    f.close()




print ('done')
