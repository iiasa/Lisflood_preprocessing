"""
# -------------------------------------------------------------------------
# Name:        Use MERIT coordinates of upstream area to create shapefiles
# Purpose:     uses upstream area of MERIT (UPA) and GRDC station data
#              to create shapefiles  from rivernetwork MERIT data
#
# Author:      PB
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022

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

# ----------------------------------------------------------------------


# import pyflwdir, some dependencies
import geopandas as gp
import numpy as np
import rasterio
from rasterio import features
import pyflwdir

import sys
import os
import os.path
import warnings

# local convenience methods (see utils.py script in notebooks folder)
#from utils import vectorize  # convenience method to vectorize rasters
#from utils import quickplot, colors, cm  # data specific quick plot method


# convenience method for vectorizing a raster
def vectorize(data, nodata, transform, name="value"):
    feats_gen = features.shapes(
        data,
        mask=data != nodata,
        transform=transform,
        connectivity=8,
    )
    feats = [
        {"geometry": geom, "properties": {name: val}} for geom, val in list(feats_gen)
    ]

    # parse to geopandas for plotting / writing to file
    gdf = gp.GeoDataFrame.from_features(feats, crs=crs)
    gdf[name] = gdf[name].astype(data.dtype)
    return gdf

#---------------------------------------------------


#-----------------------------------------------

def getBasinInd(grdc_shape,newshape):

    p_inter = gp.overlay(grdc_shape, newshape, how='intersection')
    p_union = gp.overlay(grdc_shape, newshape, how='union')

    pint_area = p_inter.area.sum()
    puni_area = p_union.area.sum()
    indicator = pint_area / puni_area

    # s1 = grdc_shape.symmetric_difference(p2, align=False)
    # sym_min = s1.area

    return indicator

# --------------------------------------------

warnings.filterwarnings("ignore")

rootshape = "../shape_glofas_3sec/"
lddname = "../data/danube_fd.tif"



sizenoarea = 5
#----------------------------------------------------------

glofas_stations = "glofas_Merit_2corr.txt"
f = open(glofas_stations, "r")
glofas = f.readlines()
f.close()
header = glofas[0]
glofas = glofas [1:]

glofas_Merit = "glofas_Merit_3.txt"
f = open(glofas_Merit, "w")
f.write(header)
f.close()

lddfile = "P:/watmodel/CWATM/cwatm_input_danube_1min/processing/routing/ldd/danube_fd.tif"
print ("load ldd")

src = rasterio.open(lddfile, "r")
flwdir = src.read(1)
transform = src.transform
latlon = src.crs.to_epsg()
crs = src.crs
src.close()

print ("make ldd")
flw = pyflwdir.from_array(flwdir, ftype='d8', transform=transform,check_ftype = False, latlon=True)
flwdir = 0

print ("done ldd")
"""
z = flw.upstream_area('km2')
name4 = "acc_danube_3sec.tif"
with rasterio.open(name4, 'w', driver='GTiff',
                   height=z.shape[0], width=z.shape[1], count=1, dtype=z.dtype,
                   crs=crs1, transform=transform1, ) as dst:
    dst.write(z, 1)

"""

for stationNo in range(0,len(glofas)):
    # ................................

    noerror = True
    
    station = glofas[stationNo].split("\t")

    upsreal = float(station[6])
    # lon,lat  -> changed because of flw library
    coord = [ float(station[21]), float(station[20])]
    # lat lon


    glofas_no =(f"{int(station[1]):04d}")
    stname = station[7] + "_"+ station[8] + "_"+ station[3]
    stname = stname.replace("ö", "oe")
    stname = stname.replace("ü", "ue")
    stname = stname.replace("ä", "ae")
    stname = stname.replace("ß", "ss")
    stname = stname.replace("Ö", "Oe")
    stname = stname.replace("Ü", "Ue")
    stname = stname.replace("Ä", "Ae")
    stname = stname.replace(")", "_")
    stname = stname.replace("(", "_")
    stname = stname.replace(" ", "_")
    provider = station[1]

    print ("-----------------------")
    print (stationNo,glofas_no,stname)


    noerror = True
    if noerror:
        print ("basin")

        try:
            subbasins = flw.basins(xy=(coord))

        except Exception as e:
            noerror = False
            print("Conversion to basin not working")

        if noerror:
            print ("vectorize")
            # vectorize subbasins using the vectorize convenience method from utils.py
            gdf_bas = vectorize(subbasins.astype(np.int32), 0, flw.transform, name=stname)

            gdf_bas['river'] = station[7]
            #gdf_bas['station'] = station[5]
            gdf_bas['country'] = station[3]
            gdf_bas['lat']= coord[1]
            gdf_bas['lon']= coord[0]
            gdf_bas['area'] = station[22]

            #gdf_bas['prov'] = station[2]
            gdf_bas['provider_no'] = station[2]
            gdf_bas['area_provider'] = station[6]
            #gdf_bas['altitude'] = station[7]
            gdf_bas['lat_provider'] = station[4]
            gdf_bas['lon_provider'] = station[5]
            #gdf_bas['startyear'] = station[10]
            #gdf_bas['endyear'] = station[11]

            shapefile1 = rootshape + glofas_no + "_" + stname + "_3sec.shp"
            shapefile1 = shapefile1.replace("ö", "oe")
            shapefile1 = shapefile1.replace("ü", "ue")
            shapefile1 = shapefile1.replace("ä", "ae")
            shapefile1 = shapefile1.replace("ß", "ss")
            shapefile1 = shapefile1.replace("Ö", "Oe")
            shapefile1 = shapefile1.replace("Ü", "Ue")
            shapefile1 = shapefile1.replace("Ä", "Ae")
            shapefile1 = shapefile1.replace("(", "_")
            shapefile1 = shapefile1.replace(")", "_")
            shapefile1 = shapefile1.replace(" ", "_")


            gdf_bas.to_file(shapefile1)
            #gdf_bas.head()
            #areashape = 0.000001 * np.sum(gdf_bas.area)


            if noerror:
                s = glofas[stationNo][:-1]
                #s = s + "\t" + str(areashape) + "\n"

                print(s)
                f = open(glofas_Merit, "a")
                f.write(s)
                f.close()
                #print ("--------- -------------------------")

    #g1 = gdf_bas.simplify(0.001)
    #g1.to_file("out_simpl21.shp")


print ("Done")