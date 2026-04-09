# Correction of station location for EFAS and GloFAS calibration

To compare simulated and observed discharge, grid-based hydrological models must fit reported station locations to the resolution-dependent gridded river network. 
In most cases this is done by comparing reported basin area for the station versus the upstream area calculated from the river network (LDD). For EFAS and GLOFAS manual correction of station also plays an important factor.

We use an intersection-over-union ratio approach to select station locations on a coarser grid-scale, reducing the errors in assigning stations to the correct upstream basin. The approaqch is explained in the method part and in Burek and Smilovic (2023).

Datasets to run the programs and the results will be placed on the ECMWF shared diskspace.

For the Danube **GloFAS** we looked at 46 stations and corrected 24 stations automatically, but only 6 with relevant changes.

For the Danube **EFAS** we looked at 315 stations and corrected 84 stations, but only 9 with relevant changes.

The small number of major corrections reflects the already high quality for station selection processes due to several rounds of EFA/Glofas calibrations.
