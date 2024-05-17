 # -*- coding: utf-8 -*-
"""Please refer to quick_guide.pdf for usage instructions"""

import os
import sys
import datetime
import time
import glob
import string

from osgeo import gdal
from osgeo import osr
from osgeo import gdalconst

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from matplotlib import rcParams
from matplotlib import rc
import pandas
from configparser import ConfigParser
import hydroStats

import warnings
warnings.filterwarnings("ignore")

#rc('font', **{'family': 'serif', 'serif': ['Palatino']})
#rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=False)
rcParams.update({'figure.autolayout': True})

t = time.time()

########################################################################
#   Read settings file
########################################################################

#iniFile = os.path.normpath(sys.argv[1])
iniFile = "settings1.txt"
#stationID = "G1005"

(drive, path) = os.path.splitdrive(iniFile)
(path, fil) = os.path.split(path)
print (">> Reading settings file (" + fil + ")...")

if not (os.path.isfile(iniFile)):
    print("No inifile found or error reading")
    sys.exit()

parser = ConfigParser()
parser.read(iniFile)

root = parser.get('DEFAULT', 'Root')
rootbasin = parser.get('DEFAULT','Rootbasin')
root = os.path.join(root, rootbasin)

#file_CatchmentsToProcess = os.path.normpath(sys.argv[2])
#file_CatchmentsToProcess = os.path.join(root, parser.get('MultiComputer', 'CatchmentsToProcess'))

ForcingStart = datetime.datetime.strptime(parser.get('DEFAULT', 'ForcingStart'), "%d/%m/%Y")  # Start of forcing
ForcingEnd = datetime.datetime.strptime(parser.get('DEFAULT', 'ForcingEnd'), "%d/%m/%Y")  # Start of forcing

#WarmupDays = int(parser.get('DEFAULT', 'WarmupDays'))
WarmupDays = 0

#CatchmentDataPath = os.path.join(root, parser.get('Path', 'CatchmentDataPath'))
SubCatchmentPath = os.path.join(root, parser.get('Path', 'SubCatchmentPath'))

#path_temp = os.path.join(root, parser.get('Path', 'Temp'))
#path_maps = os.path.join(root, os.path.join(parser.get('Path', 'CatchmentDataPath'), "maps_pcraster"))
#path_result = os.path.join(root, parser.get('Path', 'Result'))

#Qtss_csv = os.path.join(root, parser.get('ObservedData', 'Qtss'))
Qgis_csv = os.path.join(root, parser.get('ObservedData', 'Qgis'))
Qgis_out = os.path.join(root, parser.get('ObservedData', 'QgisOut'))

########################################################################
#   Loop through catchments and perform calibration
########################################################################

print (">> Reading Qgis2.csv file...")

stationdata = pandas.read_csv(os.path.join(Qgis_csv), sep=",", index_col=0)
stationdata_sorted = stationdata.sort_values(by=['CatchmentArea'], ascending=True)

#CatchmentsToProcess = pandas.read_csv(file_CatchmentsToProcess, sep=",", header=None)

# Make figures directory
figures_path = os.path.join(SubCatchmentPath, "FIGURES")
#if not os.path.exists(figures_path):
#    os.mkdir(os.path.join(SubCatchmentPath, "FIGURES"))

numberdone = 0

# Delete contents of figures directory
#for filename in glob.glob(os.path.join(SubCatchmentPath, "FIGURES", '*.*')):
#    os.remove(filename)


#Series = CatchmentsToProcess.loc[:, 0]
for index, row in stationdata.iterrows():
   #if index == int(stationID[1:]):
   #if len(Series[Series == str(row["ID"])]) == 0:  # Only process catchments whose ID is in the CatchmentsToProcess.txt file
   #    continue
   if row["ID"] == "G0137":
      iiii =1
   figfile = os.path.join(SubCatchmentPath, "FIGURES", row['ID'] + "_budyko.png")
   #while True:
   figfile = os.path.join(SubCatchmentPath, "FIGURES",row['ID']+"_budyko.png")
   if not(os.path.exists(figfile)):
     path_subcatch = os.path.join(SubCatchmentPath, row['ID'])

     maskfile = path_subcatch + "/maps/mask.map"
     nf2 = gdal.Open(maskfile, gdalconst.GA_ReadOnly)
     band = nf2.GetRasterBand(1)
     mapnp = band.ReadAsArray(0, 0, nf2.RasterXSize, nf2.RasterYSize).astype(np.float32)
     mapnp[mapnp > 1] = 0
     sum = int(np.nansum(mapnp))
     print("Number of cells: ", sum)




     ########################################################################
     #   FIGURES  FIGURES  FIGURES   FIGURES  FIGURES  FIGURES  FIGURES     #
     ########################################################################

     # FIGURE of Budyko

     filename = os.path.join(path_subcatch, "out", "09_best", "totalET_totaltot.txt")
     if os.path.exists(filename):

      print(row['ID'])  # For debugging
      filename = os.path.join(path_subcatch, "out", "00_000", "Precipitation_totaltot.txt")
      P = pandas.read_csv(filename, header=None, delim_whitespace=True, skiprows=3)
      filename = os.path.join(path_subcatch, "out", "00_000", "ETRef_totaltot.txt")
      ETP = pandas.read_csv(filename, header=None, delim_whitespace=True, skiprows=3)
      X = ETP / P

      # simulation 0
      filename = os.path.join(path_subcatch, "out", "00_000", "totalET_totaltot.txt")
      Q1 = pandas.read_csv(filename, header=None, delim_whitespace=True, skiprows=3)
      y0 = Q1 / P

      #filename = os.path.join(path_subcatch, "out", "00_002", "totalET_totaltot.txt")
      #Q1 = pandas.read_csv(filename, header=None, delim_whitespace=True, skiprows=3)
      #yNS = Q1 / P

      filename = os.path.join(path_subcatch, "out", "09_best", "totalET_totaltot.txt")
      Q1 = pandas.read_csv(filename, header=None, delim_whitespace=True, skiprows=3)
      ybest = Q1 / P

      t = np.arange(2)
      #x = np.arange(0, 2.51, 0.01)
      x = np.arange(0, 3.01, 0.01)
      maxy = x.copy()
      maxy[maxy > 1] = 1.0
      budyko = 1 + x - (1 + x ** 2.6) ** (1 / 2.6)
      budyko1 = 1 + x - (1 + x ** 2.3) ** (1 / 2.3)
      budyko2 = 1 + x - (1 + x ** 2.9) ** (1 / 2.9)



      fig = plt.figure()
      gs = plt.GridSpec(13, 6)
      plots_grids = {
          'title': (0, slice(None, None)),  # first row - Title
          'calibration': (slice(1, 5), slice(None, None)),  # second row - plot (a) calibration period
          'validation': (slice(5, 9), slice(None, None)),  # third row - plot (b) validation period
          'climatology_cal': (slice(9, 13), slice(0, 3)),  # fourth row - plot (c) monthly discharge (calibration)
          'climatology_val': (slice(9, 13), slice(3, 6)),  # fourth row - plot (d) monthly discharge (validation)
      }

      # ticks config
      months = np.array([9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # water year


      # TEXT OF CALIBRATION RESULTS
      ax0 = plt.subplot(gs[plots_grids['title']])  # first row
      ax0.set_axis_off()
      #\texts = '{}: {} at Otta'.format(row['ID'], row['Station'])
      km2 = r'$km^2$'
      texts = '{} '.format( row['ID'])
      #texts = '{} at {}, basin area: {} {}'.format( row['River'],row['Station'],str(row['DrainArPro']),km2)
      #texts_filtered = filter(lambda x: x in string.printable, texts)
      texts_filtered = texts
      ax0.text(0.5, 0.0, texts_filtered, verticalalignment='top', horizontalalignment='center', transform=ax0.transAxes, fontsize=15)

      # FIGURE OF CALIBRATION PERIOD TIME SERIES
      #Dates_Cal = Q.loc[Cal_Start:Cal_End].index
      #Q_sim_Cal = Q.loc[Cal_Start:Cal_End].iloc[:, 0].values
      #Q_obs_Cal = Q.loc[Cal_Start:Cal_End].iloc[:, 1].values
      ax1 = plt.subplot(gs[plots_grids['calibration']])  # second row
      ###max_y_cal = max(Q_sim_Cal.max(), Q_obs_Cal.max()) * 1.3
      ###ax1.set_ybound(upper=max_y_cal)
      ax1.axis('auto')
      ax1.set_adjustable('datalim')
      # format the ticks
      ###ax1.xaxis.set_major_locator(mdates.YearLocator())
      ###ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
      ###ax1.xaxis.set_minor_locator(mdates.MonthLocator())

      #ax1.plot(Dates_Cal.to_pydatetime(), Q_sim_Cal, 'r', Dates_Cal.to_pydatetime(), Q_obs_Cal, 'b')
      ax1.plot(x, budyko, 'b', linewidth=2.0)
      ax1.plot(x, maxy, 'r', linewidth=2.0)

      ax1.scatter(X, y0, s=10, facecolors='none', edgecolors='g')
      ax1.scatter(X, ybest, s=20, facecolors='none', edgecolors='r')


      #ax1.plot(x, budyko, 'b')
      #ax1.set_title('Budyko plot')
      plt.xlim([0, 3.0])
      plt.ylim([0, 1.05])
      plt.ylabel(r'E/P')
      plt.xlabel(r'$E_{pot}/P$')

      leg = ax1.legend(['Budyko curve', 'Budyko limit', 'Sim0', 'Sim Budyko'], loc=2, fancybox=True, framealpha=0.8,
                       prop={'size': 12}, labelspacing=0.1)

      leg.get_frame().set_edgecolor('white')
      plt.text(0.15, 0.40, 'demand limit', rotation=35, color='red', fontsize=10)
      plt.text(1.4, 1.04, 'supply limit', color='red', fontsize=10)
      ax1.text(0.8, 0.2,
               r'Budyco curve: $\frac{E}{P}=1+\frac{E_{pot}}{P}-(1+(\frac{E_{pot}}{P})^\omega)^{1/\omega}$  with: $\omega=2.6$',
               fontsize=12)
      # Hide the right and top spines
      ax1.spines['right'].set_visible(False)
      ax1.spines['top'].set_visible(False)
      ax1.yaxis.set_ticks_position('left')
      ax1.xaxis.set_ticks_position('bottom')

      #ax1.set_ylabel(r'$Discharge [m^3 / s]$')
      filename = os.path.join(path_subcatch, "out", "00_000", "Precipitation_totaltot.txt")
      Prec1 = pandas.read_table(filename, sep=r"\s+", index_col=0, skiprows=3, header=None, skipinitialspace=True)
      Prec = Prec1.index.values

      filename = os.path.join(path_subcatch, "out", "00_000", "ETRef_totaltot.txt")
      ETP1 = pandas.read_table(filename, sep=r"\s+", index_col=0, skiprows=3, header=None, skipinitialspace=True)
      ETP = ETP1.index.values

      filename = os.path.join(path_subcatch, "out", "09_best", "totalET_totaltot.txt")
      sim1 = pandas.read_table(filename, sep=r"\s+", index_col=0, skiprows=3, header=None, skipinitialspace=True)
      sim = sim1.index.values

      BudykoX = ETP / Prec
      BudykoY = sim / Prec
      budyko = hydroStats.Budyko(BudykoX, BudykoY)
      bb = "Number of cells: "+ str(sum) + " Budyko sum: {0:.2f}".format(budyko)
      if budyko>100:
          bb ="Number of cells: "+ str(sum) + " Budyko sum: >100"

      ax1.text(0.8, 0.1,bb  , fontsize=12)
      # budyko = HydroStats.budw1d(BudykoX, BudykoY)
      print("Budyko Sum: " + "{0:.3f}".format(budyko))
      stationdata.loc[index, 'no_cells'] = sum
      stationdata.loc[index, 'Budyko_sum'] = budyko


      # FIGURES OF CALIBRATION EVOLUTION
      #front_history = pandas.read_csv(os.path.join(path_subcatch, "front_history.csv"), sep=",", parse_dates=True, index_col=0)
      adjustprops = dict(left=0.1, bottom=0, right=1, top=0.5, wspace=-0.2, hspace=0.0)
      fig.subplots_adjust(**adjustprops)

      plt.draw()
      fig.set_size_inches(22 / 2.54, 30 / 2.54)

      figname = os.path.join(SubCatchmentPath, "FIGURES", row['ID'] + "_budyko.png")

      fig.savefig(figname, dpi=300, format='PNG', bbox_inches='tight')
      plt.close("all")
      numberdone = numberdone + 1
      ii = 1

stationdata.to_csv(Qgis_out, ',')
print ("Number done",numberdone)