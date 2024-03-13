
"""
# -------------------------------------------------------------------------
# Name:        Shapefile in low-res
# Purpose:     creates shapefiles and a list in lower resolution (5 or 30 arcmin)
#
# Author:      PB
#
# Created:     15/05/2022
# Copyright:   (c) PB 2022

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
indups: scoring on upstream area fitting between low-res area and high-res area
ups1:   low-res upstream area of best upstream fitting
indshape: scoring on similarity between low-res shape and high res shape
shapeloc: location (index of 0-24) of best fitting shape similarity
ind:    scoring result for best upstream fitting based on upstream area fitting and shape similarity
indups: scoring on upstream area fitting between low-res area and high-res area
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

import geopandas as gp
import numpy as np
import rasterio
from rasterio import features
import pyflwdir

import sys
import os
import os.path
import warnings
#-----------------------------------------------



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
#----------------------------------------------

warnings.filterwarnings("ignore")


# INPUT DATA
rootshape = "../shape_glofas_3sec/"
maketxt = "glofas_Merit_2corr.txt"

# input 5 arcmin
rootshape2 = "../shape_3min/"
rootshape3 =  "../shape_glofas/"

lddname = "ldd_3min.tif"
upsname = "upArea_repaired.nc"
# output
grdc_Merit = "basins_3min.txt"



f = open(maketxt, "r")

makeshp = f.readlines()
f.close()
header = makeshp[0].rstrip()
makeshp = makeshp[1:]

header += "\tlat_3min\tlon_3min\tarea_3min\tdifferentfromGLOFAS\n"
f = open(grdc_Merit, "w")
f.write(header)
f.close()

# read upstream area of coarse grid
src = rasterio.open(upsname)
ups = src.read(1)

# read ldd and store attributes
src = rasterio.open(lddname, "r")
flwdir = src.read(1)
transform1 = src.transform
latlon = src.crs.to_epsg()
crs = src.crs
src.close()
flwdir[flwdir<1]=5

# create river network
flw = pyflwdir.from_array(flwdir, ftype="ldd", transform=transform1, check_ftype=False, latlon=True)


#---------------------------------------

# range of 5x5 array -> this is where in the coarse grid the best point can be found

rangexy = np.arange(-2./20.,2.9/20.,1./20.0)
threshold = 100
reso ="3min"

# +++++++++++++++++++++++++++++++++
# run through als stations
notsame = 0
for stationNo in range(len(makeshp)):

    station = makeshp[stationNo].split("\t")
    upsreal = float(station[6])
    upsmodel = float(station[22])

    lat_lisf = float(station[11])
    lon_lisf = float(station[10])

    """
    try:
        alti = float(station[333])
    except:
        alti = -999.
    """

    startd = station[12]
    endd = station[13]

    coord = [ float(station[20]), float(station[21])]


    efas_no =(f"{int(station[1]):04d}")
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

    shapefile1 = rootshape + efas_no + "_" + stname + "_3sec.shp"
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

    if upsreal< 0:
        upsreal = float(station[25])

    # only station bigger equal than 1000km2 (5min) 9000km2 (30min)
    if (upsreal >= threshold) or (upsmodel >= threshold):

        if os.path.isfile(shapefile1):
            s2 = gp.GeoDataFrame.from_file(shapefile1)

            j = 0
            ind1 = []
            ind2 = []
            ind = []
            ups25 =[]
            for y in rangexy:
                for x in rangexy:
                    #print (j)
                    xx = coord[1] + x
                    yy = coord[0] + y
                    if (xx>-9999) and (xx<99999999999999999):
                        c = [xx,yy]
                        subbasins = flw.basins(xy=(c))
                        basin1km = vectorize(subbasins.astype(np.int32), 0, flw.transform, name="basin")
                        #shapefile3 = root + "5arcmin/shapetest/"+ grdc_no + "_" + str(j) + ".shp"
                        #basin1km.to_file(shapefile3)

                        # calculate union and intersection of shapes
                        p_inter = gp.overlay(s2, basin1km, how='intersection')
                        p_union = gp.overlay(s2, basin1km, how='union')
                        pint_area = p_inter.area.sum()
                        puni_area = p_union.area.sum()
                        indshape = (pint_area / puni_area)

                        # 25 values of shape similarity in ind1
                        ind1.append(indshape)

                        # get upstream area of coarse grid
                        col = int((xx - 8.0) * 20 + 0.0000001)
                        row = int((51.0 - yy) * 20 + 0.0000001)
                        upsvalue =  ups[row,col] / 1000000.

                        if upsreal == 0 or upsvalue == 0:
                            indups = 0
                        else:
                            if upsreal < upsvalue:
                                indups = upsreal / upsvalue
                            else:
                                indups = upsvalue / upsreal

                        # 25 values of area difference in ind2
                        ind2.append(indups)
                        # 25 values of area ups25
                        ups25.append(upsvalue)

                        # calculate lenght to origin (0,0)
                        # multiobjective ups and shape
                        #ind.append(np.sqrt((1-indups)**2+(1-indshape)**2))

                        j += 1

            #---------------------------------
            # maximum of shape simi and ups accordance
            maxups = np.max(ind2)
            maxshape = np.max(ind1)
            # index 0-24 -> of max shape and max ups
            shapemax = np.where(ind1 == maxshape)[0][0]

            upsshape = ups25[shapemax]
            ups12 = ups25[12]
            diffups = abs(upsshape-ups12)
            
            # if difference <=50 km2 than use middle point
            if diffups <= 50.0:
                 pzt= 100 * abs(1 - ups12/upsshape)
                 # only if percentage of change is equal or smaller than 1 pct
                 if pzt <= 1.0:
                    #print(shapemax, upsshape, ups12, abs(upsshape - ups12))
                    maxshape = ind1[12]
                    shapemax = 12


            upsmax = np.where(ind2 == maxups)[0][0]

            # using multiobjective ups and shape
            #minind = np.min(ind)
            #indmin = np.where(ind==minind)[0][0]

            # using only shape
            indmin = shapemax

            # get location in a 5x5 matrix
            y= indmin // 5
            x = indmin % 5

            # resulting coordinate , but still a multiply of the 3sec coordinate
            yy = coord[0] + rangexy[y]
            xx = coord[1] + rangexy[x]

            # getting the loc and the upstream area on coarse resolution
            col = int((xx - 8.0) * 20  + 0.0000001)
            row = int((51.0 - yy) * 20 + 0.0000001)
            upsvalue = ups[row, col] / 1000000
            locx =   col / 20.0 + 8.0 + 1./40. #xx  //60.0 * 60.0+1/120
            locy = 51.0 - row / 20 - 1./40.    #yy  #//60. * 60.0+ 1/120

            col_lisf = int((lon_lisf - 8.0) * 20 + 0.0000001)
            row_lisf = int((51.0 - lat_lisf) * 20 + 0.0000001)
            locxlisf = col_lisf / 20.0 + 8.0 + 1./40.
            locylisf = 51.0 - row_lisf / 20 - 1. / 40.

            nots = "0"
            if (locxlisf != locx) or (locylisf != locy):
                print ("Lisflood loc different",notsame)
                notsame = notsame + 1
                nots = "1"

                subbasins = flw.basins(xy=([lon_lisf,lat_lisf]))
                gdf_bas = vectorize(subbasins.astype(np.int32), 0, flw.transform, name=stname)
                shapefile3 = rootshape3 + efas_no  + "_glofas_" + stname + "_3min.shp"

                gdf_bas['river'] = station[7]
                #gdf_bas['station'] = station[5]
                gdf_bas['country'] = station[3]
                gdf_bas['no_glofas'] = station[1]
                gdf_bas['provider_no'] = station[2]
                gdf_bas['area_provider'] = station[6]
                #gdf_bas['altitude'] = station[7]
                gdf_bas['lat_provider'] = station[4]
                gdf_bas['lon_provider'] = station[5]

                gdf_bas['area_new'] = upsvalue
                gdf_bas['lat_new'] = locy
                gdf_bas['lon_new'] = locx

                gdf_bas['area_glofas'] = station[9]
                gdf_bas['lat_glofas']= lat_lisf
                gdf_bas['lon_glofas']= lon_lisf

                gdf_bas.to_file(shapefile3)



            subbasins = flw.basins(xy=([xx,yy]))
            gdf_bas = vectorize(subbasins.astype(np.int32), 0, flw.transform, name=stname)
            shapefile2 = rootshape2 + efas_no  + "_" + stname + "_3min.shp"

            gdf_bas['river'] = station[7]
            #gdf_bas['station'] = station[5]
            gdf_bas['country'] = station[3]
            gdf_bas['no_GLOFAS'] = station[1]
            gdf_bas['provider_no'] = station[2]
            gdf_bas['area_provider'] = station[6]
            #gdf_bas['altitude'] = station[7]
            gdf_bas['lat_provider'] = station[4]
            gdf_bas['lon_provider'] = station[5]

            gdf_bas['lat_3sec']= coord[1]
            gdf_bas['lon_3sec']= coord[0]
            gdf_bas['area_3sec'] = station[22]

            gdf_bas['area_3min'] = upsvalue
            gdf_bas['lat'] = locy
            gdf_bas['lon'] = locx

            gdf_bas.to_file(shapefile2)


            s =  makeshp[stationNo][:-1]
            s = s + "\t" + f"{locy:8.5f}" + "\t" + f"{locx:6.5f}" + "\t" + f"{upsvalue:9.0f}"


            # add if Glofas is different
            s = s + "\t" + nots
            print (s)
            s = s + "\n"
            f = open(grdc_Merit, "a")
            f.write(s)
            f.close()

print ("done")
print (notsame)