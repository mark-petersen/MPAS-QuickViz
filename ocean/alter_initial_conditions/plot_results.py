#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = './'
initFileName = 'init.nc'
outputFileName = 'output.nc'

############################## domain and IC
nx = 32
ny = 4
import numpy as np
dc = 1.0
Lx = nx*dc
Ly = ny*dc*np.sqrt(3)/2
kx = 1
ky = 0
nVertLevels=1

print('loading libraries...')
from datetime import date
import matplotlib.pyplot as plt
import xarray as xr
initDS = xr.open_dataset(runDir+initFileName)
outDS = xr.open_dataset(runDir+outputFileName)

nCells = initDS.dims['nCells']
xCell = initDS.variables['xCell']
yCell = initDS.variables['yCell']
nEdges = initDS.dims['nEdges']
xEdge = initDS.variables['xEdge']
yEdge = initDS.variables['yEdge']
angleEdge = initDS.variables['angleEdge']
nVertices = initDS.dims['nVertices']
xVertex = initDS.variables['xVertex']
yVertex = initDS.variables['yVertex']

#divergenceSol = initDS.variables['divergenceSol']
#relativeVorticitySol = initDS.variables['relativeVorticitySol']
#del2GradDivVelocitySol = initDS.variables['del2GradDivVelocitySol']
#del2GradVortVelocitySol = initDS.variables['del2GradVortVelocitySol']
#
#divergence = initDS.variables['divergence']
#relativeVorticity = initDS.variables['relativeVorticity']
#del2GradDivVelocityTendency = initDS.variables['del2GradDivVelocityTendency']
#del2GradVortVelocityTendency = initDS.variables['del2GradVortVelocityTendency']

k=0

figdpi = 300
fig = plt.figure(figsize=(40,12))
varNames = ['divergenceSol', 'relativeVorticitySol','del2GradDivVelocitySol','del2GradVortVelocitySol',
            'divergence', 'relativeVorticity','del2GradDivVelocityTendency','del2GradVortVelocityTendency']
loc = ['cell','vertex','edge','edge','cell','vertex','edge','edge']
f = ['in','in','in','in','out','out','out','out']

for j in range(len(varNames)):
    ax = plt.subplot(2,4,j+1)

    if f[j]=='in':
        var = initDS.variables[varNames[j]][:,k]
    elif f[j]=='out':
        var = outDS.variables[varNames[j]][0,:,k]
    if loc[j]=='cell':
        im = plt.scatter(xCell/1000,yCell/1000,c=var,s=16,marker='H',cmap=plt.cm.jet)
    elif loc[j]=='edge':
        im = plt.scatter(xEdge/1000,yEdge/1000,c=var,s=16,marker='s',cmap=plt.cm.jet)
    elif loc[j]=='vertex':
        im = plt.scatter(xVertex/1000,yVertex/1000,c=var,s=16,marker='^',cmap=plt.cm.jet)
    plt.colorbar(im, ax=ax)
    if j>=6:
        varSol = initDS.variables[varNames[j-4]][:,k]
        diff = var - varSol
        print('nx={} '.format(np.sqrt(nCells))+'maxabs: {:9.2E}'.format(np.max(abs(diff))),varNames[j])
        #print('rms: {:9.2E}'.format(np.mean(diff**2)))
        plt.title(varNames[j]+' max diff: {:9.2E}'.format(np.max(diff))+' nx={}'.format(np.sqrt(nCells)))
    else:
        plt.title(varNames[j])

    #plt.xlabel('x, km')
    #plt.ylabel('y, km')
figfile = 'plot_del2_nx{}'.format(np.sqrt(nCells))+'.png'
plt.savefig(figfile) #, bbox_inches='tight')
plt.close()

