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

############################## domain and IC
print('loading libraries...')
from datetime import date
import matplotlib.pyplot as plt
import xarray as xr

for p in range(4,10):
    N = 2**p
    Lx = 1024.0e3
    nx = N
    ny = N
    import numpy as np
    dc = int(Lx/N)
    Ly = ny*dc*np.sqrt(3)/2
    nVertLevels=1
    fileName = 'base_mesh_{}x{}.nc'.format(nx,ny)
    newFileName = 'init_{}x{}.nc'.format(nx,ny)
    #print('planar_hex --nx {} --ny {} --dc {} -o '.format(nx,ny,dc)+fileName)

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
    del2VelocitySol = np.zeros([nEdges,nVertLevels])
    
    k=0
    IC=1
    # Create initial conditions for sin(x)*sin(y)
    if IC==1:
        kux=1; kuy=2; kvx=2; kvy=3; pi2 = 2.0*np.pi
        u = np.sin( kux*pi2/Lx*xEdge[:] ) * \
            np.sin( kuy*pi2/Ly*yEdge[:] )
        v = np.sin( kvx*pi2/Lx*xEdge[:] ) * \
            np.sin( kvy*pi2/Ly*yEdge[:] )
        ux= kux*pi2/Lx* \
            np.cos( kux*pi2/Lx*xCell[:] ) * \
            np.sin( kuy*pi2/Ly*yCell[:] )
        vy= kvy*pi2/Ly* \
            np.sin( kvx*pi2/Lx*xCell[:] ) * \
            np.cos( kvy*pi2/Ly*yCell[:] )
        uy= kuy*pi2/Ly* \
            np.sin( kux*pi2/Lx*xVertex[:] ) * \
            np.cos( kuy*pi2/Ly*yVertex[:] )
        vx= kvx*pi2/Lx* \
            np.cos( kvx*pi2/Lx*xVertex[:] ) * \
            np.sin( kvy*pi2/Ly*yVertex[:] )
        uxy = kux*pi2/Lx * kuy*pi2/Ly * \
            np.cos( kux*pi2/Lx*xEdge[:] ) * \
            np.cos( kuy*pi2/Ly*yEdge[:] )
        vxy = kvx*pi2/Lx * kvy*pi2/Ly * \
            np.cos( kvx*pi2/Lx*xEdge[:] ) * \
            np.cos( kvy*pi2/Ly*yEdge[:] )
        uxx = -(kux*pi2/Lx)**2*u
        uyy = -(kuy*pi2/Ly)**2*u
        vxx = -(kvx*pi2/Lx)**2*v
        vyy = -(kvy*pi2/Ly)**2*v
    if IC==2:
        # Create initial conditions for uyy and vxx only:
        kux=0; kuy=1; kvx=2; kvy=0; pi2 = 2.0*np.pi
        u = np.sin( kuy*pi2/Ly*yEdge[:] )
        v = np.sin( kvx*pi2/Lx*xEdge[:] )
        ux= kux*pi2/Lx* \
            np.cos( kux*pi2/Lx*xCell[:] ) * \
            np.sin( kuy*pi2/Ly*yCell[:] )
        vy= kvy*pi2/Ly* \
            np.sin( kvx*pi2/Lx*xCell[:] ) * \
            np.cos( kvy*pi2/Ly*yCell[:] )
        uy= kuy*pi2/Ly* \
            np.cos( kuy*pi2/Ly*yVertex[:] )
        vx= kvx*pi2/Lx* \
            np.cos( kvx*pi2/Lx*xVertex[:] ) 
        uxy = kux*pi2/Lx * kuy*pi2/Ly * \
            np.cos( kux*pi2/Lx*xEdge[:] ) * \
            np.cos( kuy*pi2/Ly*yEdge[:] )
        vxy = kvx*pi2/Lx * kvy*pi2/Ly * \
            np.cos( kvx*pi2/Lx*xEdge[:] ) * \
            np.cos( kvy*pi2/Ly*yEdge[:] )
        uxx = -(kux*pi2/Lx)**2*u
        uyy = -(kuy*pi2/Ly)**2*u
        vxx = -(kvx*pi2/Lx)**2*v
        vyy = -(kvy*pi2/Ly)**2*v
    zonalVelocityEdge[:,k] = u
    meridionalVelocityEdge[:,k] = v
    normalVelocity[:,k] = \
          np.cos(angleEdge[:]) * u \
        + np.sin(angleEdge[:]) * v
    
    # Create exact solutions:
    
    divergenceSol[:,k] = ux + vy
    relativeVorticitySol[:,k] = vx - uy
    del2GradDivVelocitySol[:,k] = \
          np.cos(angleEdge[:]) * (uxx + vxy) \
        + np.sin(angleEdge[:]) * (uxy + vyy)
    del2GradVortVelocitySol[:,k] = \
          np.cos(angleEdge[:]) * (uyy - vxy) \
        + np.sin(angleEdge[:]) * (vxx - uxy)
    del2VelocitySol[:,k] = \
          np.cos(angleEdge[:]) * (uxx + uyy) \
        + np.sin(angleEdge[:]) * (vxx + vyy)
    
    print('ln -isf '+newFileName+' init.nc; srun -n 1 ocean_model; python plot_results.py >> results.txt')
    #mesh.expand_dims({'nVertLevels':nVertLevels})
    mesh["zonalVelocityEdge"]=(['nEdges','nVertLevels'], zonalVelocityEdge)
    mesh["meridionalVelocityEdge"]=(['nEdges','nVertLevels'], meridionalVelocityEdge)
    mesh["normalVelocity"]=(['nEdges','nVertLevels'], normalVelocity)
    mesh["divergenceSol"]=(['nCells','nVertLevels'], divergenceSol)
    mesh["relativeVorticitySol"]=(['nVertices','nVertLevels'], relativeVorticitySol)
    mesh["del2GradDivVelocitySol"]=(['nEdges','nVertLevels'], del2GradDivVelocitySol)
    mesh["del2GradVortVelocitySol"]=(['nEdges','nVertLevels'], del2GradVortVelocitySol)
    mesh["del2VelocitySol"]=(['nEdges','nVertLevels'], del2VelocitySol)
    
    H = 1000.0
    maxLevelCell = np.ones([nCells],dtype=np.int32); mesh["maxLevelCell"]=(['nCells'], maxLevelCell)
    refBottomDepth = H*np.ones([nVertLevels]); mesh["refBottomDepth"]=(['nVertLevels'], refBottomDepth)
    refZMid = H/2*np.ones([nVertLevels]); mesh["refZMid"]=(['nVertLevels'], refZMid)
    layerThickness = H*np.ones([nCells,nVertLevels]); mesh["layerThickness"]=(['nCells','nVertLevels'], layerThickness)
    restingThickness = H*np.ones([nCells,nVertLevels]); mesh["restingThickness"]=(['nCells','nVertLevels'], restingThickness)
    vertCoordMovementWeights = np.ones([nVertLevels],dtype=np.int32); mesh["vertCoordMovementWeights"]=(['nVertLevels'], vertCoordMovementWeights)
    
    mesh.to_netcdf(path=runDir+newFileName)

exit()

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
