#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
"""
Plot vertical profiles from individual columns
Mark Petersen
September 2022
"""

############################## model files, run dirs
runDir = '/lustre/scratch5/turquoise/mpeterse/runs/s08a/'
runDir = './'
fileName = 'base_mesh.nc'
newFileName = 'init.nc'

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
print('planar_hex --nx {} --ny {} --dc {} -o base_mesh_{}x{}.nc'.format(nx,ny,dc,nx,ny))

print('loading libraries...')
from datetime import date
import matplotlib.pyplot as plt
import xarray as xr
#from netCDF4 import Dataset
mesh = xr.open_dataset(runDir+fileName)

nCells = mesh.dims['nCells']
xCell = mesh.variables['xCell']
yCell = mesh.variables['yCell']
nEdges = mesh.dims['nEdges']
xEdge = mesh.variables['xEdge']
yEdge = mesh.variables['yEdge']
angleEdge = mesh.variables['angleEdge']
nVertices = mesh.dims['nVertices']
xVertex = mesh.variables['xVertex']
yVertex = mesh.variables['yVertex']

zonalVelocityEdge = np.zeros([nEdges,nVertLevels])
meridionalVelocityEdge = np.zeros([nEdges,nVertLevels])
normalVelocity = np.zeros([nEdges,nVertLevels])
divergenceSol = np.zeros([nCells,nVertLevels])
relativeVorticitySol = np.zeros([nVertices,nVertLevels])
del2GradDivVelocitySol = np.zeros([nEdges,nVertLevels])
del2GradVortVelocitySol = np.zeros([nEdges,nVertLevels])

k=0
zonalVelocityEdge[:,k] = np.sin( 2*np.pi*kx / Lx * xEdge[:] )
meridionalVelocityEdge[:,k] = np.sin( 2*np.pi*ky / Ly * yEdge[:] )
normalVelocity[:,k] = np.cos(angleEdge[:]) * zonalVelocityEdge[:,k] + np.sin(angleEdge[:]) * meridionalVelocityEdge[:,k]

divergenceSol[:,k] = 2*np.pi*kx / Lx * np.cos( 2*np.pi*kx / Lx * xCell[:] )
relativeVorticitySol[:,k] = 0.0 * xVertex[:]
del2GradDivVelocitySol[:,k] = -(2*np.pi*kx / Lx)**2 * np.sin( 2*np.pi*kx / Lx * xEdge[:] )
del2GradVortVelocitySol[:,k] = 0.0 * xEdge[:]

print('write file:')
#mesh.expand_dims({'nVertLevels':nVertLevels})
mesh["zonalVelocityEdge"]=(['nEdges','nVertLevels'], zonalVelocityEdge)
mesh["meridionalVelocityEdge"]=(['nEdges','nVertLevels'], meridionalVelocityEdge)
mesh["normalVelocity"]=(['nEdges','nVertLevels'], normalVelocity)
mesh["divergenceSol"]=(['nCells','nVertLevels'], divergenceSol)
mesh["relativeVorticitySol"]=(['nVertices','nVertLevels'], relativeVorticitySol)
mesh["del2GradDivVelocitySol"]=(['nEdges','nVertLevels'], del2GradDivVelocitySol)
mesh["del2GradVortVelocitySol"]=(['nEdges','nVertLevels'], del2GradVortVelocitySol)

H = 1000.0
maxLevelCell = np.ones([nCells],dtype=np.int32); mesh["maxLevelCell"]=(['nCells'], maxLevelCell)
refBottomDepth = H*np.ones([nVertLevels]); mesh["refBottomDepth"]=(['nVertLevels'], refBottomDepth)
refZMid = H/2*np.ones([nVertLevels]); mesh["refZMid"]=(['nVertLevels'], refZMid)
layerThickness = H*np.ones([nCells,nVertLevels]); mesh["layerThickness"]=(['nCells','nVertLevels'], layerThickness)
restingThickness = H*np.ones([nCells,nVertLevels]); mesh["restingThickness"]=(['nCells','nVertLevels'], restingThickness)
vertCoordMovementWeights = np.ones([nVertLevels],dtype=np.int32); mesh["vertCoordMovementWeights"]=(['nVertLevels'], vertCoordMovementWeights)

divergenceSol = np.zeros([nCells,nVertLevels])
relativeVorticitySol = np.zeros([nVertices,nVertLevels])
del2GradDivVelocitySol = np.zeros([nEdges,nVertLevels])
del2GradVortVelocitySol = np.zeros([nEdges,nVertLevels])
mesh.to_netcdf(path=runDir+newFileName)

figdpi = 300
fig = plt.figure(figsize=(20,12))
varNames = ['zonalVelocityEdge', 'meridionalVelocityEdge','normalVelocity','angleEdge']
for j in range(len(varNames)):
    ax = plt.subplot(2,2,j+1)
    var = mesh.variables[varNames[j]]
    im = plt.scatter(xEdge/1000,yEdge/1000,c=var,s=15,marker='s',cmap=plt.cm.jet)
    plt.scatter(xCell/1000,yCell/1000,c='k',s=8,marker='H')
    #plt.gca().invert_yaxis()
    #plt.grid()
    plt.set_cmap('jet')
    plt.colorbar(im, ax=ax)
    plt.title(varNames[j])
    plt.xlabel('x, km')
    plt.ylabel('y, km')
figfile = 'plot_init.png'
plt.savefig(figfile) #, bbox_inches='tight')
plt.close()
