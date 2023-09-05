#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s06o/ocean/global_ocean/QU30/PHC/init/initial_state/'
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s07a/ocean/global_ocean/WC14/PHC/init/initial_state/'
fileName = 'initial_state.nc'
newFileName = 'initial_state_sloped.nc'

#import matplotlib as mpl
from datetime import date
import numpy as np
#import matplotlib.pyplot as plt
import xarray as xr
from netCDF4 import Dataset

deg2rad = 3.14159/180.0
rad2deg = 180.0/3.14159


print('read: '+runDir+fileName)
mesh = xr.open_dataset(runDir+fileName)
nCells = mesh.dims['nCells']
latCell = mesh.variables['latCell']
lonCell = mesh.variables['lonCell']
maxLevelCell = mesh.variables['maxLevelCell']
bottomDepth = mesh.variables['bottomDepth']
refBottomDepth = mesh.variables['refBottomDepth']

#cellList = np.where(np.logical_and(
#    latCell>latMin*deg2rad, latCell<=(latMin+latWidth)*deg2rad))[0]

# Southern strip
latWidth = 2
latMin = 9
latMax = latMin + latWidth
lonMin = -63 + 360
lonMax =  360
cellList = np.where(np.logical_and(np.logical_and(np.logical_and(
    latCell>latMin*deg2rad, latCell<=latMax*deg2rad), 
    lonCell>lonMin*deg2rad), lonCell<=lonMax*deg2rad))[0]
print('before iCell loop south')
for iCell in cellList:
    maxLevelCell[iCell] = min(maxLevelCell[iCell], 
        5 + (60-5)*(latCell[iCell]*rad2deg - latMin)/(latMax - latMin) )
    bottomDepth[iCell] = min(bottomDepth[iCell], 
        refBottomDepth[maxLevelCell[iCell]-1])

# Northern strip
latWidth = 2
latMax = 43
latMin = latMax - latWidth
lonMin = -80 + 360
lonMax =  360
cellList = np.where(np.logical_and(np.logical_and(np.logical_and(
    latCell>latMin*deg2rad, latCell<=latMax*deg2rad), 
    lonCell>lonMin*deg2rad), lonCell<=lonMax*deg2rad))[0]
print('before iCell loop north')
for iCell in cellList:
    maxLevelCell[iCell] = min(maxLevelCell[iCell], 
        60 - (60-5)*(latCell[iCell]*rad2deg - latMin)/(latMax - latMin) )
    bottomDepth[iCell] = min(bottomDepth[iCell], 
        refBottomDepth[maxLevelCell[iCell]-1])

#mesh["cullCell"]=(['nCells'],  cullCell)
mesh.to_netcdf(path=runDir+newFileName)
