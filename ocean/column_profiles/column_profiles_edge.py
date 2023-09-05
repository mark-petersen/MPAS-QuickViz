#!/usr/bin/env python
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s/'
simName = 's03n'
fileName = '/output.nc'
meshName = '/init.nc'

lonMin = 0.193
lonMax = 0.210
latMin = 0.477
latMax = 0.490

#from __future__ import absolute_import, division, print_function, \
#    unicode_literals
#import os
#import glob
import matplotlib as mpl
#mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as cols
#from matplotlib.pyplot import cm
#from matplotlib.colors import BoundaryNorm
#import cmocean
import xarray as xr
from netCDF4 import Dataset
#from mpas_analysis.shared.io.utility import decode_strings
#import gsw

print('reading data from:')
print(runDir+simName+meshName)
print(runDir+simName+fileName)
mesh = xr.open_dataset(runDir+simName+meshName)
data = xr.open_dataset(runDir+simName+fileName)

z = mesh.refBottomDepth.values
nVertLevels = len(z)
nEdges = data.dims['nEdges']

latEdge = mesh.variables['latEdge']
lonEdge = mesh.variables['lonEdge']
maxLevelEdge = data.variables['maxLevelEdgeTop']

edgeList = np.where(np.logical_and(np.logical_and(np.logical_and(latEdge>latMin, latEdge<latMax), lonEdge>lonMin), lonEdge<lonMax))[0]

iTime = 0
fig = plt.figure(figsize=(20,12))
varNames = ['divergence','vertVelocityTop','vertTransportVelocityTop','temperature', 'salinity','density','pressure','zMid']
varNames = ['temperature', 'salinity']
varNames = ['divergence','vertVelocityTop','temperature', 'salinity']
varNames = ['normalVelocity']
for j in range(len(varNames)):
    plt.subplot(2,2,j+1)
    var = data.variables[varNames[j]]
    for i in range(len(edgeList)):
        iEdge = int(edgeList[i])
        k = int(maxLevelEdge[iEdge])
        #varData = var[iTime,iEdge,0:k]
        varData = np.log10(abs(var[iTime,iEdge,0:k]))
        #varData = np.where(varData>-1e20,varData,np.NAN)
        plt.plot(varData,np.arange(1,k+1),label='edge '+str(iEdge))
    plt.gca().invert_yaxis()
    plt.title(varNames[j])
    plt.grid()
    plt.legend()

figfile = 'vert_profiles_' +simName+ '.png'
plt.savefig(figfile, bbox_inches='tight')
plt.close()
