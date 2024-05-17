#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calibration tool for Hydrological models
using a distributed evolutionary algorithms in python
DEAP library
https://github.com/DEAP/deap/blob/master/README.md

Félix-Antoine Fortin, François-Michel De Rainville, Marc-André Gardner, Marc Parizeau and Christian Gagné, "DEAP: Evolutionary Algorithms Made Easy", Journal of Machine Learning Research, vol. 13, pp. 2171-2175, jul 2012

The calibration tool was created by Hylke Beck 2014 (JRC, Princeton) hylkeb@princeton.edu
Thanks Hylke for making it available for use and modification
Modified by Peter Burek

The submodule Hydrostats was created 2011 by:
Sat Kumar Tomer (modified by Hylke Beck)
Please see his book "Python in Hydrology"   http://greenteapress.com/pythonhydro/pythonhydro.pdf

"""
import os
import sys
import shutil
import hydroStats
import array
import random
import numpy as np
import datetime
import glob
from deap import algorithms
from deap import base
from deap import benchmarks
from deap import creator
from deap import tools
import pandas

from osgeo import gdal
from osgeo import osr
from osgeo import gdalconst

import multiprocessing
import time
from configparser import ConfigParser
import glob
from subprocess import Popen, PIPE

import ast
from sys import platform
import pickle
import platform as plat

from multiprocessing import Process, freeze_support
#if sys.version_info[1] > 7: #test if python > 3.7
from multiprocessing import shared_memory

#import importlib
import run_cwatm
from copy import copy, deepcopy

import warnings
import signal

## Set global parameter
global gen
gen = 0
WarmupDays = 0

#--------------------------------------------
def signal_handler(sig, frame):

	# This function will handle the interrupt signal (Ctrl+C)
	print(" =============== Exiting gracefully...")
	try:
		cname = os.getenv('COMPUTERNAME')
		shm1 = shared_memory.SharedMemory(name=cname)
		print("--- open 1")
	except:
		print("no open 1")
	try:
		shm1.close()
		print("--- close 1")
	except:
		print("no close 1")
	try:
		shm1.unlink()
		print("--- unlink 1")
	except:
		print("no unlink 1")
	# ---------------------------
	# --------------------------------------

	try:
		shm5.close()
		rint("--- close 5")
	except:
		print("no close 5")
	try:
		shm5.unlink()
		print("--- unlink 5")
	except:
		print("no unlink 5")
	try:
		shm6.close()
		print("--- close 6")
	except:
		print("no close 6")
	try:
		shm6.unlink()
		print("--- unlink 6")
	except:
		print("no unlink 6")

	sys.exit(0)
# ---------------------------------
# FIRST RUN FOR Meteo

def firstmeteorun():
	id = 0
	run_rand_id = str(id // 1000).zfill(2) + "_" + str(id % 1000).zfill(3)
	directory_run = os.path.join(path_subcatch)
	gaugeloc = str(station['lon']) + " " + str(station['lat'])
	template_xml_new = template_xml
	template_xml_new = template_xml_new.replace("%SubCatchmentPath", path_subcatch)
	#template_xml_new = template_xml_new.replace("%meteoData", path_meteoData)
	template_xml_new = template_xml_new.replace("%root", root)

	template_xml_new = template_xml_new.replace('%gaugeloc', gaugeloc)  # Gauge location
	#template_xml_new = template_xml_new.replace('%CalStart', Cal_Realstart)  # Date of Cal starting
	#template_xml_new = template_xml_new.replace('%CalSpin', Cal_Start1)  # Date of Cal starting
	#template_xml_new = template_xml_new.replace('%CalEnd', Cal_End1)

	#template_xml_new = template_xml_new.replace('%inflowflag', inflowflag)
	#template_xml_new = template_xml_new.replace('%inflowDir', path_inflow)
	#template_xml_new = template_xml_new.replace('%inflowpoints', inflowloc)
	#template_xml_new = template_xml_new.replace('%inflowtss', "inflow.tss")
	
	#template_xml_new = template_xml_new.replace('%initDanube', "Danube_1km_19990110.nc")

	for ii in range(0, len(ParamRanges) - 1):
		template_xml_new = template_xml_new.replace("%" + ParamRanges.index[ii], "1.0")
	# replace output directory
	template_xml_new = template_xml_new.replace('%run_rand_id', directory_run)
	settings = os.path.join(directory_run, os.path.basename(ModelSettings_template)[:-4] + '-meteo.ini')
	f = open(settings,"w")
	f.write(template_xml_new)
	f.close()
	return settings

# -----------------------------------
def writetime_user(listPC,stationID,gen):
	# write time and user

	h = plat.uname()[1]
	t = datetime.datetime.now()
	tt = t.strftime("%Y-%m-%d %H:%M")
	s = str(tt) + "," + h +  "," + str(gen) + "\n"

	file = os.path.join(listPC, stationID + "_list.txt")
	if os.path.isfile(file):
		f = open(file, 'a')
	else:
		f = open(file, 'w')
	f.write(s)
	f.close()


########################################################################
#   Read settings file
########################################################################

#iniFile = os.path.normpath(sys.argv[1])
#stationID = os.path.normpath(sys.argv[2])

iniFile = "P:/watmodel/calibration/global5min_budyko/settings1.txt"
stationID = "G1005"
shared = True
#print ("calibration_single2c")
#print (iniFile, stationID)

sstart =2010
ssend = 2019
sspin = 3

# work around because .exe do not take the right sys.argv after first time
# works on python 3.8 (shared memory is new): dtype <U68 -> improve at the moment the settings file name has to be the same number of chars
# freeze_support is important for .exe translation of multiprocessing
if shared:
	try:
		cname = os.getenv('COMPUTERNAME')
		shm1 =shared_memory.SharedMemory(name=cname)
		saveshared = False
	except:
		saveshared = True


	if saveshared:
	#if os.path.isfile(iniFile):
		info = np.array([iniFile.ljust(180)+","+stationID.ljust(19)])
		cname = os.getenv('COMPUTERNAME')
		shm1 = shared_memory.SharedMemory(name=cname,create=True, size=info.nbytes)
		# Now create a NumPy array backed by shared memory
		b1 = np.ndarray(info.shape, dtype=info.dtype, buffer=shm1.buf)
		b1[:] = info[:]  # Copy the original data into shared memory


	else:
		ii =1
		cname = os.getenv('COMPUTERNAME')
		i1 = np.array(["bla".ljust(180) + "," + "bla".ljust(19)])
		shm1 =shared_memory.SharedMemory(name=cname)
		c1 = np.ndarray((2), dtype=i1.dtype, buffer=shm1.buf)
		c1 = c1[0].split(",")
		iniFile = c1[0].rstrip()
		stationID = c1[1].rstrip()

if not(os.path.isfile(iniFile)):
	print(iniFile)
	print("No inifile found or error reading")
	sys.exit()

#print (iniFile, stationID)

parser = ConfigParser()
parser.read(iniFile)

if platform == "win32":
	root = parser.get('DEFAULT','Root')
else:
	root = parser.get('DEFAULT','RootLinux')

rootbasin = parser.get('DEFAULT','Rootbasin')
root = os.path.join(root, rootbasin)
#print (root)

try:
	ForcingStart = datetime.datetime.strptime(parser.get('DEFAULT', 'ForcingStart'),"%d/%m/%Y %H:%M")  # Start of forcing
	ForcingEnd = datetime.datetime.strptime(parser.get('DEFAULT', 'ForcingEnd'), "%d/%m/%Y %H:%M")  # Start of forcing
except:
	ForcingStart = datetime.datetime.strptime(parser.get('DEFAULT', 'ForcingStart'),"%d/%m/%Y")  # Start of forcing
	ForcingEnd = datetime.datetime.strptime(parser.get('DEFAULT', 'ForcingEnd'), "%d/%m/%Y")

timeperiod = parser.get('DEFAULT','timeperiod')
if timeperiod == "monthly":
	monthly = 1
	dischargetss = 'totalET_totaltot.txt'
	frequen = 'MS'
else:
	monthly = 0
	dischargetss = 'totalET_totaltot.txt'
	frequen = 'd'

Qtss_csv = os.path.join(root,parser.get('ObservedData', 'Qtss'))
Qgis_csv = os.path.join(root,parser.get('ObservedData', 'Qgis'))

path_result = os.path.join(root,parser.get('Path', 'Result'))
ParamRangesPath = os.path.join(root,parser.get('Path','ParamRanges'))
SubCatchmentPath = os.path.join(root,parser.get('Path','SubCatchmentPath'))
#path_meteoData = parser.get('Path','MeteoData')

#Qtss_csv = os.path.join(rootbasin,parser.get('ObservedData', 'Qtss'))
#Qtss_col = parser.get('ObservedData', 'Column')

modeltemplate = os.path.join(root,parser.get('Path','Templates'))
ModelSettings_template = os.path.join(modeltemplate,parser.get('Templates','ModelSettings'))
RunModel_template = os.path.join(root,parser.get('Templates','RunModel'))

spinoffyears = int(parser.get('DEFAULT','SpinoffYears'))


# Multi computer as executable runs
if platform == "win32":
	Run = parser.get('MultiComputer', 'RunCwatm')
	Run1 = Run.split(" ")
	if len(Run1) > 1:
		RunCwatm = Run1[0] + " " + os.path.join(root, Run1[1])
	else:
		RunCwatm = os.path.join(root, Run1[0])
	set = ModelSettings_template.split(".")
	ModelSettings_template = set[0] +".ini"


else:
	Run1 = parser.get('MultiComputer', 'RunCwatmLinux')
	RunCwatm = Run1.split(" ")[0] + " " + os.path.join(root,Run1.split(" ")[1])
	set = ModelSettings_template.split(".")
	ModelSettings_template = set[0] +"Linux.ini"



listPC = os.path.join(root,parser.get('MultiComputer', 'listPC'))


use_multiprocessing = int(parser.get('DEAP','use_multiprocessing'))

try:
	pool_limit = int(parser.get('DEAP','pool_limit'))
except:
	pool_limit = 10000

ngen = int(parser.get('DEAP','ngen'))
mu = int(parser.get('DEAP','mu'))
lambda_ = int(parser.get('DEAP','lambda_'))
maximize =  parser.getboolean('DEAP','maximize')
if maximize: maxDeap = 1.0
else: maxDeap = -1.0

try:
	select_best = int(parser.get('DEAP','select_best'))
except:
	select_best = lambda_
try:
	start_from_gen1 = parser.getboolean('DEAP','start_from_gen1')
except:
	start_from_gen1 = False

firstrun = parser.getboolean('Option', 'firstrun')   # using default run as first run
bestrun = parser.getboolean('Option', 'bestrun')

########################################################################
#   Preparation for calibration
########################################################################

# Load xml and .bat template files
runmodel = os.path.splitext(RunModel_template)[0]



f = open(ModelSettings_template,"r")
template_xml = f.read()
f.close()

# Load parameter range file
ParamRanges = pandas.read_csv(ParamRangesPath,sep=",",index_col=0)

stationdata = pandas.read_csv(Qgis_csv, sep=",", index_col=0)
station = stationdata.loc[stationdata["ID"]==stationID].iloc[0]

Cal_S = datetime.datetime.strptime("01/01/"+str(sstart), '%d/%m/%Y')
#Cal_S = datetime.datetime.strptime(station['Cal_Start'], '%d/%m/%Y')
Cal_Start = Cal_S.strftime('%Y-%m-%d')
Cal_Start1 = Cal_S.strftime('%d/%m/%Y')

Cal_End = datetime.datetime.strptime("31/12/"+str(ssend), '%d/%m/%Y').strftime('%Y-%m-%d')
Cal_End1 = datetime.datetime.strptime("31/12/"+str(ssend), '%d/%m/%Y').strftime('%d/%m/%Y')

#Cal_End = datetime.datetime.strptime(station['Cal_End'], '%d/%m/%Y').strftime('%Y-%m-%d')
#Cal_End1 = datetime.datetime.strptime(station['Cal_End'], '%d/%m/%Y').strftime('%d/%m/%Y')
Cal_Realstart = datetime.datetime(Cal_S.year - sspin, Cal_S.month, Cal_S.day, 0, 0).strftime('%d/%m/%Y')

path_subcatch = os.path.join(SubCatchmentPath, station['ID'])

# Load observed streamflow
"""
streamflow_data = pandas.read_csv(Qtss_csv, sep=",", parse_dates=True,dayfirst = True, index_col=0)
observed_streamflow = streamflow_data[station['ID']]
#observed_streamflow = observed_streamflow[ForcingStart:ForcingEnd]
observed_streamflow = observed_streamflow[Cal_Start:Cal_End]
observed_streamflow[observed_streamflow<-900]= np.nan

path_inflow = os.path.join(SubCatchmentPath, station['ID'],"inflow")

if station["Inflow"] > 0:
	inflowflag = "True"
	inflowloc_txt = os.path.join(path_inflow, "inflowloc2.txt")
	f = open(inflowloc_txt, "r")
	inflowloc = f.read()
	f.close()

else:
	inflowflag = "False"
	inflowloc = "0 0"

"""


# first standard parameter set
# Snowmelt, crop KC, soil depth,pref. flow, arno beta, groundwater recession, runoff conc., routing, manning, No of run
# recalculated to a population setting
if firstrun:
	para_first2 = []
	for ii in range(0, len(ParamRanges) - 1):
		delta = float(ParamRanges.iloc[ii, 1]) - float(ParamRanges.iloc[ii, 0])
		if delta == 0:
			para_first2.append(0.)
		else:
			para_first2.append((float(ParamRanges.iloc[ii, 2]) - float(ParamRanges.iloc[ii, 0])) / delta)

ii = 1

# Check last run and repeat:
maskfile = path_subcatch + "/maps/mask.map"
nf2 = gdal.Open(maskfile, gdalconst.GA_ReadOnly)
band = nf2.GetRasterBand(1)
mapnp = band.ReadAsArray(0, 0, nf2.RasterXSize, nf2.RasterYSize).astype(np.float32)
mapnp[mapnp>1] = 0
sum = np.nansum(mapnp)
print ("Number of cells: ", sum)


writetime_user(listPC, stationID,0)

bestrunnotexists = True
bestrundir = path_subcatch +"/out/"+ str(ngen+1).zfill(2)+"_best"
if os.path.isdir(bestrundir):
	#for f in os.listdir(bestrundir):
	#	os.remove(os.path.join(bestrundir, f))
	#time.sleep(1)
	bestrunnotexists = False

if bestrunnotexists:
	if shared:
		if saveshared:

			fsettingsfile = firstmeteorun()
			#meteo = np.array([1,2,3])
			meteo,success, last_dis = run_cwatm.main(fsettingsfile, ['-lk'])
			meteo = meteo.astype('float32')
			meteoshape = np.array([meteo.shape])
			shm5 = shared_memory.SharedMemory(name=stationID+'meteoshape', create=True, size=meteoshape.nbytes)
			b2 = np.ndarray(meteoshape.shape, dtype=meteoshape.dtype, buffer=shm5.buf)
			b2[:] = meteoshape[:]

			shm6 = shared_memory.SharedMemory(name=stationID+'meteodata', create=True, size=meteo.nbytes)
			b3 = np.ndarray(meteo.shape, dtype=meteo.dtype, buffer=shm6.buf)
			b3[:] = meteo[:]

		else:
			shm5 =shared_memory.SharedMemory(name=stationID+'meteoshape')
			c2 = np.ndarray((3), dtype='int32', buffer=shm5.buf)

			shm6 =shared_memory.SharedMemory(name=stationID+'meteodata')
			meteo = np.ndarray((c2[0],c2[1],c2[2]), dtype='float32', buffer=shm6.buf)
	else:
		meteo = []

#	file = os.path.join(listPC, stationID + "_list.txt")
#os.remove(file)
#print ("Something not working on ",stationID)






########################################################################
#   Function for running the model, returns objective function scores
########################################################################

def RunModel(Individual):

	# Convert scaled parameter values ranging from 0 to 1 to usncaled parameter values
	Parameters = [None] * len(ParamRanges)
	for ii in range(0,len(ParamRanges-1)):
		Parameters[ii] = Individual[ii]*(float(ParamRanges.iloc[ii,1])-float(ParamRanges.iloc[ii,0]))+float(ParamRanges.iloc[ii,0])

	# Note: The following code must be identical to the code near the end where the model is run
	# using the "best" parameter set. This code:
	# 1) Modifies the settings file containing the unscaled parameter values amongst other things
	# 2) Makes a .bat file to run the model
	# 3) Run the model and loads the simulated streamflow

	# Random number is appended to settings and .bat files to avoid simultaneous editing
	#run_rand_id = str(gen).zfill(2) + "_" + str(int(random.random()*100000000)).zfill(10)
	id =int(Individual[-1])
	run_rand_id = str(id//1000).zfill(2) + "_" + str(id%1000).zfill(3)

	ii = id%1000
	if ii // 20 == ii / 20:
		writetime_user(listPC, stationID,str(id//1000).zfill(2))

	directory_run = os.path.join(path_subcatch,"out", run_rand_id)

	gaugeloc = str(station['lon']) + " " + str(station['lat'])

	template_xml_new = template_xml

	template_xml_new = template_xml_new.replace("%SubCatchmentPath", path_subcatch)
	#template_xml_new = template_xml_new.replace("%meteoData", path_meteoData)
	template_xml_new = template_xml_new.replace("%root", root)

	template_xml_new = template_xml_new.replace('%gaugeloc', gaugeloc)  # Gauge location
	#template_xml_new = template_xml_new.replace('%CalStart', Cal_Realstart)  # Date of Cal starting
	#template_xml_new = template_xml_new.replace('%CalSpin', Cal_Start1)  # Date of Cal starting
	#template_xml_new = template_xml_new.replace('%CalEnd', Cal_End1)

	#template_xml_new = template_xml_new.replace('%inflowflag', inflowflag)
	#template_xml_new = template_xml_new.replace('%inflowDir', path_inflow)
	#template_xml_new = template_xml_new.replace('%inflowpoints', inflowloc)
	#template_xml_new = template_xml_new.replace('%inflowtss', "inflow.tss")

	#template_xml_new = template_xml_new.replace('%initDanube', "Danube_1km_19991231.nc")

	for ii in range(0,len(ParamRanges)-1):
		template_xml_new = template_xml_new.replace("%"+ParamRanges.index[ii],str(Parameters[ii]))
	# replace output directory
	template_xml_new = template_xml_new.replace('%run_rand_id', directory_run)

	if os.path.isdir(directory_run):
		if os.path.exists(os.path.join(directory_run,dischargetss)):
			runmodel = False
		else:
			runmodel = True
			shutil.rmtree(directory_run)
	else: runmodel = True


	if runmodel:
		os.mkdir(directory_run)
		settings = os.path.join(directory_run,os.path.basename(ModelSettings_template)[:-4] + '-Run' + run_rand_id + '.ini')
		f = open(settings, "w")
		f.write(template_xml_new)
		f.close()

		#template_bat_new = template_bat
		# python P:/watmodel/CWATM\cwatm_input_1km_pinzgau\CWATM_160421_lucaModflow/run_cwatm.py %run -v
		if use_multiprocessing == 0:
			template_bat_new = RunCwatm + " %run -l\npause"
		else:
			template_bat_new = RunCwatm + " %run -v"
		template_bat_new = template_bat_new.replace('%run',os.path.basename(ModelSettings_template)[:-4]+'-Run'+run_rand_id+'.ini')
		#runfile = os.path.join(directory_run,os.path.basename(RunModel_template)[:-4]+run_rand_id) + ".bat"

		if platform == "win32":
			runfile = os.path.join(directory_run, os.path.basename(RunModel_template)[:-4] + run_rand_id) + ".bat"
		else:
			runfile = os.path.join(directory_run, os.path.basename(RunModel_template)[:-4] + run_rand_id) + ".sh"

		f = open(runfile, "w")
		f.write(template_bat_new)
		f.close()

		if shared:
			if use_multiprocessing == 1:
				success, last_dis = run_cwatm.mainwarm(settings, ['-v'], meteo)
			else:
				success, last_dis = run_cwatm.mainwarm(settings, ['-l'], meteo)
		else:
			if use_multiprocessing == 1:
				success, last_dis = run_cwatm.mainwarm(settings, ['-v'], meteo)
			else:
				success, last_dis = run_cwatm.mainwarm(settings, ['-l'], meteo)

	#if not(shared):
	# P = observed_data/Precipitation_totaltot.txt
	# ETP = observed_data/ETRef_totaltot.txt
	p_csv = os.path.join(directory_run, "Precipitation_totaltot.txt")
	etp_csv = os.path.join(directory_run, "ETRef_totaltot.txt")

	Prec1 = pandas.read_table(p_csv, sep=r"\s+", index_col=0, skiprows=3, header=None, skipinitialspace=True)
	ETP1 = pandas.read_table(etp_csv, sep=r"\s+", index_col=0, skiprows=3, header=None, skipinitialspace=True)
	Prec = Prec1.index.values
	ETP = ETP1.index.values
	BudykoX = ETP / Prec


	Qsim_tss = os.path.join(directory_run,dischargetss)

	if os.path.isfile(Qsim_tss)==False:
		print("run_rand_id: "+str(run_rand_id)+" File: "+ Qsim_tss)
		raise Exception("No actual evaporation found. Probably the model failed to start? Check the log files of the run!")
	try:
		#simulated_streamflow = pandas.read_csv(Qsim_tss,sep=r"\s+",index_col=0,skiprows=4,header=None,skipinitialspace=True)
		# Compute objective function score
		sim1 = pandas.read_table(Qsim_tss, sep=r"\s+", index_col=0, skiprows=3, header=None, skipinitialspace=True)
		sim = sim1.index.values


	except:
		shutil.rmtree(directory_run)
		raise Exception("Can't read simulated streamflow. It is deleted and run again!")


	BudykoY = sim / Prec
	budyko = hydroStats.Budyko(BudykoX, BudykoY)
	#budyko = HydroStats.budw1d(BudykoX, BudykoY)
	print ("   run_rand_id: "+str(run_rand_id)+", Budyko Sum: "+"{0:.3f}".format(budyko))

	s = str(run_rand_id) + ", " + "{0:.4f}".format(budyko) + "\n"

	runlogfile = os.path.join(path_subcatch, "runs_log.csv")

	if os.path.isfile(runlogfile):
		with open(runlogfile, 'r') as file:
			lines = file.readlines()

		found = False
		for i, line in enumerate(lines):
			if str(run_rand_id) in line:
				lines[i] = s
				found = True
				break

		if found:
			with open(runlogfile, 'w') as file:
				file.writelines(lines)
		else:
			with open(runlogfile, 'a') as file:
				file.write(s)
	else:
		with open(runlogfile, 'w') as file:
			file.write("run_no, Budyko\n")
			file.write(s)

	return budyko, # If using just one objective function, put a comma at the end!!!

	"""
	COR = hydroStats.correlation(s=Qsim,o=Qobs,warmup=WarmupDays)
	print("   run_rand_id: "+str(run_rand_id)+", COR "+"{0:.3f}".format(COR))
	with open(os.path.join(path_subcatch,"runs_log.csv"), "a") as myfile:
		myfile.write(str(run_rand_id)+","+str(COR)+"\n")
	return COR, # If using just one objective function, put a comma at the end!!!


	NSE = hydroStats.NS(s=Qsim, o=Qobs, warmup=WarmupDays)
	print "   run_rand_id: " + str(run_rand_id) + ", NSE: " + "{0:.3f}".format(NSE)
	with open(os.path.join(path_subcatch, "runs_log.csv"), "a") as myfile:
		myfile.write(str(run_rand_id) + "," + str(NSE) + "\n")
	return NSE,  # If using just one objective function, put a comma at the end!!!
	"""

########################################################################
#   Perform calibration using the DEAP module
########################################################################

creator.create("FitnessMin", base.Fitness, weights=(maxDeap,))
#creator.create("Individual", array.array, typecode='d', fitness=creator.FitnessMin)
creator.create("Individual", array.array, typecode='d', fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_float", random.uniform, 0, 1)

# Structure initializers
toolbox.register("Individual", tools.initRepeat, creator.Individual, toolbox.attr_float, len(ParamRanges))
toolbox.register("population", tools.initRepeat, list, toolbox.Individual)


def checkBounds(min, max):
	def decorator(func):
		def wrappper(*args, **kargs):
			offspring = func(*args, **kargs)
			for child in offspring:
				for i in range(len(child)):
					if child[i] > max:
						child[i] = max
					elif child[i] < min:
						child[i] = min
			return offspring
		return wrappper
	return decorator

toolbox.register("evaluate", RunModel)
toolbox.register("mate", tools.cxBlend, alpha=0.15)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.3, indpb=0.3)
toolbox.register("select", tools.selNSGA2)





history = tools.History()

toolbox.decorate("mate", checkBounds(0, 1))
toolbox.decorate("mutate", checkBounds(0, 1))

if __name__ == "__main__":

	gen = 0
	#print ("blablabla")
	freeze_support()
	t = time.time()
	writetime_user(listPC, stationID, gen)

	if bestrunnotexists==False: use_multiprocessing = False
	# if bestrun exist then no multiprocesss

	if use_multiprocessing==True:
		pool_size = multiprocessing.cpu_count() * 1
		print(pool_size, pool_limit)
		if pool_size > pool_limit: pool_size = pool_limit
		pool = multiprocessing.Pool(processes=pool_size)
		toolbox.register("map", pool.map)
		print(pool_size)


	# For someone reason, if sum of cxpb and mutpb is not one, a lot less Pareto optimal solutions are produced
	cxpb = 0.7  # The probability of mating two individuals
	mutpb = 0.3 # The probability of mutating an individual.

	effavg = np.zeros((ngen+1))*np.NaN

	if use_multiprocessing == 0:
		print ("Start calibration")


	startlater = False
	checkpoint = os.path.join(path_subcatch ,"checkpoint.pkl")
	if os.path.exists(os.path.join(checkpoint)):
		with open(checkpoint, "rb" ) as cp_file:
			cp = pickle.load(cp_file)
			populationall = cp["populationall"]
			population = cp["population"]
			start_gen = cp["generation"]
			random.setstate(cp["rndstate"])
			if start_gen > 0:
				offspring = cp["offspring"]
				halloffame =  cp["halloffame"]
				startlater = True
				gen = start_gen
		cp_file.close()

		if start_from_gen1:
			population = populationall[0].copy()
			populationall = {}
			populationall[0] = population.copy()
			startlater = True
			gen = 1
			halloffame = tools.ParetoFront()
			halloffame.update(population)
			population[:] = toolbox.select(population, select_best)



	# if population is not saved before
	else:
		population = toolbox.population(n=mu)
		# Numbering of runs
		for ii in range(mu):
			population[ii][-1]= float(gen * 1000 + ii+1)

		#first run parameter set:
		if firstrun:
			for ii in range(len(population[0])-1):
				population[0][ii] = para_first2[ii]
			population[0][-1] = 0.

	if not (startlater):
		halloffame = tools.ParetoFront()
		history.update(population)

		# saving population
		populationall ={}
		populationall[0] = population.copy()
		cp = dict(populationall=populationall, population=population, generation=gen, rndstate=random.getstate())
		with open(checkpoint, "wb") as cp_file:
			pickle.dump(cp, cp_file)
		cp_file.close()


		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in population if not ind.fitness.valid]
		fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)

		writetime_user(listPC, stationID,gen)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit

		halloffame.update(population)

		# delete all Modflow files for each run in this generation, because they are so large
		fileList = glob.glob(path_subcatch + "/out/**/Bhima*", recursive=True)
		for filePath in fileList:
			try:
				os.remove(filePath)
			except:
				print("Error while deleting file")
		time.sleep(6)


		# Loop through the different objective functions and calculate some statistics
		# from the Pareto optimal population


		effavg[0] = halloffame[0].fitness.values[0]

		gen = 0
		print(">> gen: "+str(gen)+", KGE: "+"{0:.3f}".format(effavg[gen]))
		#history.update(population)
		population[:] = toolbox.select(population, lambda_)

		# select best  population - child number (lambda_) from initial population
		gen = 1







	# Begin the generational process
	# from gen 1 .....
	conditions = {"ngen" : False, "StallFit" : False}
	while not any(conditions.values()):
		if startlater == False:
			# Vary the population
			offspring = algorithms.varOr(population, toolbox, lambda_, cxpb, mutpb)
			#offspring = algorithms.varAnd(population, toolbox, cxpb, mutpb)

			# put in the number of run
			for ii in range(lambda_):
				offspring[ii][-1] = float(gen * 1000 + ii + 1)


		# saving population
		populationall[gen] = population.copy()
		cp = dict(populationall=populationall,population=population, generation=gen, rndstate=random.getstate(), offspring=offspring, halloffame=halloffame)
		with open(checkpoint, "wb") as cp_file:
			pickle.dump(cp, cp_file)
		cp_file.close()
		startlater = False


		# ---------------------------------
		# write time and user
		writetime_user(listPC, stationID,gen)

		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit
		"""
		# delete alll Modflow files for each run in this generation, because they are so large
		fileList = glob.glob(path_subcatch + "/out/**/Bhima*", recursive=True)
		for filePath in fileList:
			try:
				os.remove(filePath)
			except:
				print("Error while deleting file")
		"""



		# Update the hall of fame with the generated individuals
		if halloffame is not None:
			halloffame = tools.ParetoFront()
			halloffame.update(offspring)

		# Select the next generation population
		population[:] = toolbox.select(population + offspring,select_best)
		history.update(population)

		# put in the number of run
		#for ii in xrange(mu):
		#	population[ii][-1] = float(gen * 1000 + ii + 1)

		# Loop through the different objective functions and calculate some statistics
		# from the Pareto optimal population
		effavg[gen] = halloffame[0].fitness.values[0]

		print(">> gen: "+str(gen)+", KGE: "+"{0:.3f}".format(effavg[gen]))

		# Terminate the optimization after ngen generations
		if gen >= ngen:
			print(">> Termination criterion ngen fulfilled.")
			conditions["ngen"] = True

		gen += 1
		# Copied and modified from algorithms.py eaMuPlusLambda until here




	# Finito
	if use_multiprocessing == True:
		pool.close()
	#elapsed = time.time() - t
	#print(">> Time elapsed: "+"{0:.2f}".format(elapsed)+" s")


	########################################################################
	#   Save calibration results
	########################################################################
	"""
	import networkx as nx

	print ("make history graph")
	graph = nx.DiGraph(history.genealogy_tree)
	graph = graph.reverse()
	# Make the graph top-down
	color =[]
	for i in graph:

		color.append(toolbox.evaluate(history.genealogy_history[i])[0])

	colors = [toolbox.evaluate(history.genealogy_history[i])[0] for i in graph]
	nx.draw(graph, node_color=color)
	graph.add_node(1, color='red', style='filled', fillcolor='blue', shape='square')

	A = nx.nx_agraph.to_agraph(graph)

	A.node_attr['style'] = 'filled'
	A.node_attr['shape'] = 'circle'
	A.node_attr['gradientangle'] = 90


	for i in A.nodes():
		n = A.get_node(i)
		n.attr['fillcolor'] = 'green;0.5:yellow'
		n.attr['fillcolor'] = color

	A.layout(prog='dot')
	A.draw('test.png',prog='dot')

	"""



	# -----------------------------------------------------------------

	# Save history of the change in objective function scores during calibration to csv file
	print(">> Saving optimization history (front_history.csv)")
	front_history = pandas.DataFrame({'gen':list(range(gen)),
									  'effavg_R':effavg[:],
									  })
	front_history.to_csv(os.path.join(path_subcatch,"front_history.csv"),',')
	# as numpy  numpy.asarray  ; numpy.savetxt("foo.csv", a, delimiter=","); a.tofile('foo.csv',sep=',',format='%10.5f')

	# Compute overall efficiency scores from the objective function scores for the
	# solutions in the Pareto optimal front
	# The overall efficiency reflects the proximity to R = 1, NSlog = 1, and B = 0 %

	front = np.array([ind.fitness.values for ind in halloffame])
	effover = 1 - np.sqrt((1-front[:,0]) ** 2)
	best = np.argmax(effover)

	# Convert the scaled parameter values of halloffame ranging from 0 to 1 to unscaled parameter values
	paramvals = np.zeros(shape=(len(halloffame),len(halloffame[0])))
	paramvals[:] = np.NaN
	for kk in range(len(halloffame)):
		for ii in range(len(ParamRanges)):
			paramvals[kk][ii] = halloffame[kk][ii]*(float(ParamRanges.iloc[ii,1])-float(ParamRanges.iloc[ii,0]))+float(ParamRanges.iloc[ii,0])

	# Save Pareto optimal solutions to csv file
	# The table is sorted by overall efficiency score
	print(">> Saving Pareto optimal solutions (pareto_front.csv)")
	ind = np.argsort(effover)[::-1]
	pareto_front = pandas.DataFrame({'effover':effover[ind],'R':front[ind,0]})
	for ii in range(len(ParamRanges)):
		pareto_front["param_"+str(ii).zfill(2)+"_"+ParamRanges.index[ii]] = paramvals[ind,ii]
	pareto_front.to_csv(os.path.join(path_subcatch,"pareto_front.csv"),',')

	# Select the "best" parameter set and run Model for the entire forcing period
	Parameters = paramvals[best,:]


	if bestrun:
		if not(os.path.exists(os.path.join(path_subcatch, "streamflow_simulated_best.tss"))):

			print(">> Running Model using the \"best\" parameter set")
			# Note: The following code must be identical to the code near the end where Model is run
			# using the "best" parameter set. This code:
			# 1) Modifies the settings file containing the unscaled parameter values amongst other things
			# 2) Makes a .bat file to run Model
			# 3) Runs Model and loads the simulated streamflow
			# Random number is appended to settings and .bat files to avoid simultaneous editing
			best_run = str(int(history.genealogy_history[1][9]))
			bl = len(best_run)

			print (history.genealogy_history[1])
			print (best_run)
			best_run = best_run[:bl - 3].zfill(2) + "_"+ best_run[bl-3:]
			directory_best = os.path.join(path_subcatch, "out", best_run)


			run_rand_id = str(gen).zfill(2) + "_best"
			template_xml_new = template_xml

			directory_run = os.path.join(path_subcatch, "out", run_rand_id)
			gaugeloc = str(station['lon']) + " " + str(station['lat'])


			#Cal_Realstart1 = datetime.datetime(Cal_S.year - spinoffyears*2, Cal_S.month, Cal_S.day, 0, 0).strftime('%d/%m/%Y')
			#Cal_Realstart1 = datetime.datetime(Cal_S.year - spinoffyears * 2, Cal_S.month, Cal_S.day, 0, 0).strftime('%d/%m/%Y')
			#Cal_Realstart_US = datetime.datetime(Cal_S.year - spinoffyears * 2, Cal_S.day, Cal_S.month, 0, 0).strftime('%d/%m/%Y')
			#Cal_Realstart1 = datetime.datetime(Cal_S.year - spinoffyears * 3, Cal_S.month, Cal_S.day, 0, 0).strftime('%d/%m/%Y')
			#Cal_Realstart_US = datetime.datetime(Cal_S.year - spinoffyears * 3, Cal_S.day, Cal_S.month, 0, 0).strftime('%d/%m/%Y')
			Cal_Realstart1 = datetime.datetime(sstart,1,1, 0, 0).strftime('%d/%m/%Y')
			Cal_Spin1 = datetime.datetime(sstart-sspin,1,1, 0, 0).strftime('%d/%m/%Y')
			Cal_Realstart_US = datetime.datetime(sstart,1,1, 0, 0).strftime('%d/%m/%Y')
			Cal_End1 = datetime.datetime(ssend,12,31, 0, 0).strftime('%d/%m/%Y')

			template_xml_new = template_xml_new.replace("%SubCatchmentPath", path_subcatch)
			#template_xml_new = template_xml_new.replace("%meteoData", path_meteoData)
			template_xml_new = template_xml_new.replace("%root", root)

			template_xml_new = template_xml_new.replace('%gaugeloc', gaugeloc)  # Gauge location
			#template_xml_new = template_xml_new.replace('%CalStart', Cal_Realstart1)  # Date of Cal starting
			#template_xml_new = template_xml_new.replace('%CalSpin', Cal_Spin1)  # Date of Cal starting
			#template_xml_new = template_xml_new.replace('%CalSpin', "None")

			#template_xml_new = template_xml_new.replace('%CalEnd', Cal_End1)

			#template_xml_new = template_xml_new.replace('%inflowflag', inflowflag)
			#template_xml_new = template_xml_new.replace('%inflowDir', path_inflow)
			#template_xml_new = template_xml_new.replace('%inflowpoints', inflowloc)
			#template_xml_new = template_xml_new.replace('%inflowtss', "inflow_last_run.tss")

			for ii in range(0,len(ParamRanges)):
				template_xml_new = template_xml_new.replace("%"+ParamRanges.index[ii],str(Parameters[ii]))
			template_xml_new = template_xml_new.replace('%run_rand_id', directory_run)

			#template_xml_new = template_xml_new.replace('%initDanube', "Danube_1km_19991231.nc")

			if bestrunnotexists:
				os.mkdir(directory_run)

			#template_xml_new = template_xml_new.replace('%InitModel',"1")
			settings = os.path.join(directory_run,os.path.basename(ModelSettings_template)[:-4]+'-Run'+run_rand_id+'.ini')
			f = open(settings, "w")
			f.write(template_xml_new)
			f.close()
			"""
			template_bat_new = template_bat
			template_bat_new = template_bat_new.replace('%run',os.path.basename(ModelSettings_template)[:-4]+'-Run'+run_rand_id+'.ini')
			"""

			if use_multiprocessing == 0:
				template_bat_new = RunCwatm + " %run -l\npause"
			else:
				template_bat_new = RunCwatm + " %run -v"
			template_bat_new = template_bat_new.replace('%run',os.path.basename(ModelSettings_template)[:-4]+'-Run'+run_rand_id+'.ini')
			runfile = os.path.join(directory_run,os.path.basename(RunModel_template)[:-4]+run_rand_id) + ".bat"

			runfile = os.path.join(directory_run, os.path.basename(RunModel_template)[:-4] + run_rand_id)
			if platform == "win32":
				runfile = runfile + ".bat"
			else:
				runfile = runfile + ".sh"


			f = open(runfile, "w")
			f.write(template_bat_new)
			f.close()

			writetime_user(listPC, stationID,gen)

			currentdir = os.getcwd()
			os.chdir(directory_run)
			"""
			p = Popen(runfile, shell=True, stdout=PIPE, stderr=PIPE, bufsize=16*1024*1024)
			output, errors = p.communicate()
			f = open("log"+run_rand_id+".txt",'w')
			content = "OUTPUT:\n"+str(output)+"\nERRORS:\n"+str(errors)
			f.write(content)
			f.close()
			os.chdir(currentdir)
			"""

			#print (directory_run)
			#print (directory_best)
			#Qsim_tss = os.path.join(directory_run, dischargetss)
			# use already existing best run
			Qsim_tss = os.path.join(directory_run, dischargetss)

			# copy files

			#shutil.copy(os.path.join(directory_best, dischargetss), os.path.join(directory_run, dischargetss))
			#shutil.copy(os.path.join(directory_best, "ETRef_totaltot.txt"), os.path.join(directory_run, "ETRef_totaltot.txt"))
			#shutil.copy(os.path.join(directory_best, "Precipitation_totaltot.txt"),os.path.join(directory_run, "Precipitation_totaltot.txt"))

			if not(os.path.exists(Qsim_tss)):
				success, last_dis = run_cwatm.mainwarm(settings, ['-l'], meteo)

			# even if monthly daily data are store here
			#dischargetss = 'discharge_daily.tss'
			#Qsim_tss = os.path.join(directory_run,dischargetss)

			simulated_streamflow = pandas.read_csv(Qsim_tss,sep=r"\s+",index_col=0,skiprows=3,header=None,skipinitialspace=True)
			#simulated_streamflow[0][simulated_streamflow[0]==1e31] = np.nan
			#Qsim = simulated_streamflow[0].values
			Qsim = simulated_streamflow.index.to_numpy()

			# Save simulated streamflow to disk
			print(">> Saving \"best\" simulated streamflow (streamflow_simulated_best.tss)")
			#Qsim = pandas.DataFrame(data=Qsim, index=pandas.date_range(ForcingStart, periods=len(Qsim), freq=frequen))
			#Qsim = pandas.DataFrame(data=Qsim, index=pandas.date_range(Cal_Realstart1, periods=len(Qsim), freq=frequen))
			Qsim = pandas.DataFrame(data=Qsim, index=pandas.date_range(Cal_Realstart_US, periods=len(Qsim), freq="d"))
			Qsim.to_csv(os.path.join(path_subcatch,"streamflow_simulated_best.tss"),',',header="")
			try: os.remove(os.path.join(path_subcatch,"out",'streamflow_simulated_best.tss'))
			except: pass
			#os.rename(Qsim_tss, os.path.join(path_subcatch,"out",'streamflow_simulated_best.tss'))

			print ("================DONE == " + stationID + " =====================================")
			print ("-------------------------------------------------------------------------------")
			print ("")

	"""
	# delete alll Modflow files for each run in this generation, because they are so large
	fileList = glob.glob(path_subcatch + "/out/**/Bhima*", recursive=True)
	for filePath in fileList:
		try:
			os.remove(filePath)
		except:
			print("Error while deleting file")
	"""

	try:
		cname = os.getenv('COMPUTERNAME')
		shm1 = shared_memory.SharedMemory(name=cname)
		print("--- open 1")
	except:
		print("no open 1")
	try:
		shm1.close()
		print("--- close 1")
	except:
		print("no close 1")
	try:
		shm1.unlink()
		print("--- unlink 1")
	except:
		print("no unlink 1")

	# ---------------------------

	try:
		shm5.close()
		rint("--- close 5")
	except:
		print("no close 5")
	try:
		shm5.unlink()
		print("--- unlink 5")
	except:
		print("no unlink 5")
	try:
		shm6.close()
		print("--- close 6")
	except:
		print("no close 6")
	try:
		shm6.unlink()
		print("--- unlink 6")
	except:
		print("no unlink 6")


	print ("delete files")

	# Delete all .xml, .bat, .tmp, and .txt files created for the runs
	for out1 in range(ngen+1):
		if out1 == 0:
			range2 = mu + 1
		else:
			range2 = lambda_ + 1
		for out2 in range(range2):
			dir1 = "%02d_" % (out1)
			dir2 = "%03d" % (out2)

			outname = os.path.join(path_subcatch, "out", dir1+dir2)
			if not((dir1+dir2)== "00_000"):
				if os.path.exists(outname):
					shutil.rmtree(outname)
					time.sleep(0.1)

	print ("Deleted all runs but best and 0")