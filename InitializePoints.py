#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:42:15 2024

@author: vedhasyamuvva
"""

# Section 2 Part 1: Distributing particles according to a multi-variate Gaussian with a chosen center, semimajor axis a, and axis ratios b/a and c/a. 

import numpy as np
import matplotlib.pyplot as plt

""" Generates N random points according to a gaussian distribution set to the chosen parameters
 
Inputs:
 ------
 center : NumPy array
     The center point of the gaussian distribution

 a : NumPy value
     The length of the semimajor axis

 ba : NumPy value
     The axis ratio of the second axis to the semimajor axis

 ca : NumPy value
     The axis ratio of the third axis to the semimajor axis
 
 N : NumPy integer
     The number of points to be generated, the default is 32x32x32 points

 Returns:
 -------
 points: Numpy array
    N random points generated with the chosen parameters

"""
def initializeGaussianPoints(center, a, ba, ca, N = 32**3):
    mean = np.array(center)
    cov = np.diag(np.array([a**2, (a*ba)**2, (a*ca)**2]))
    
    points = np.random.multivariate_normal(mean, cov, N)
    return points

""" Plots the points generated with initialize Gaussian Points
 
Inputs:
 ------
 arr : NumPy array
     array of N 3D points
"""
def plotInitialPoints(arr):
    
    ax = plt.figure().add_subplot(projection='3d')
    arr = arr.T
    ax.scatter(arr[0], arr[1], arr[2], '.', label='initial points')
    
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    plt.title('3D Scatter Plot of Gaussian Distributed Points')
    
    plt.savefig("basicpoints.png")
    plt.show()

""" Calculates the Density field for every point and returns the density field
 
Inputs:
 ------
 center : NumPy array
     The point at which the points have been centered around

 a : NumPy value
     The length of the semimajor axis

 ba : NumPy value
     The axis ratio of the second axis to the semimajor axis

 ca : NumPy value
     The axis ratio of the third axis to the semimajor axis
 
 arr : NumPy array
     the array of points previously generated randomly
     
 grid_size : NumPy value
     the number of desired points on the mesh grid
     
 Returns:
 -------
 densityField: Numpy array
    the density of of the each cell arranged in a somewhat convoluted fashion that needs to be fixed

 x: Numpy array
    array of vertices on the x-axis
   
 y: Numpy array
    array of vertices on the y-axis
    
 z: Numpy array
    array of vertices on the z-axis
"""  
def CreateDensityField(center, a, ba, ca, arr, grid_size):
    #FOR DEBUGGING PURPOSES:
    DEBUG = False
    #Create a a NxNxN Grid of points to calculate the density for 
    
    # print("center: ", center)
    # print("arr: ", arr)
    x = np.linspace(-grid_size/2, grid_size/2, grid_size+1, endpoint=True) + center[0]
    y = np.linspace(-grid_size/2, grid_size/2, grid_size+1, endpoint=True) + center[1]
    z = np.linspace(-grid_size/2, grid_size/2, grid_size+1, endpoint=True) + center[2]
    # print("x: ", x, np.shape(x))
    # print("y: ", y, np.shape(y))
    # print("z: ", z, np.shape(z))
    
    #Initialize Density Field
    
    densityField = np.zeros((grid_size+1, grid_size+1, grid_size+1))
    # print("initialized Density Field: ", densityField, np.shape(densityField))
    
    # Calculate the addition to the density field from each individual particle
    for ptcl in arr:
        #Find the grid point closest to the particle
        
        # print("particle: ", ptcl)
        xi = np.searchsorted(x, ptcl[0])
        yi = np.searchsorted(y, ptcl[1])
        zi = np.searchsorted(z, ptcl[2])
        if DEBUG: print("xi,yi,zi: ", xi,yi,zi)
        #Check to see which vertices of the box will matter
        if xi == 0: xarr = np.array([xi])
        elif xi == grid_size+1: xarr = np.array([xi-1])
        else: xarr = np.array([xi-1,xi])
        if DEBUG: print("xarr: ", xarr)
        
        if yi == 0: yarr = np.array([yi])
        elif yi >= grid_size+1: yarr = np.array([yi-1])
        else: yarr = np.array([yi-1,yi])
        if DEBUG: print("yarr: ", yarr)
        
        if zi == 0: zarr = np.array([zi])
        elif zi >= grid_size+1: zarr = np.array([zi-1])
        else: zarr = np.array([zi-1,zi])
        if DEBUG: print("zarr: ", zarr)
        
        
        #Check the vertices of the box around it
        for i in xarr:
            for j in yarr:
                for k in zarr:
                    if DEBUG: print("ijk: ", i, j, k)
                    point = np.array([x[i],y[j],z[k]])
                    if DEBUG: print("ptcl: ", ptcl)
                    if DEBUG: print("point: ", point)
                    distance = (ptcl-point)
                    if DEBUG: print("distance: ", distance)
                    dx = max(0, 1 - np.abs(distance[0]))
                    dy = max(0, 1 - np.abs(distance[1]))
                    dz = max(0, 1 - np.abs(distance[2]))
                    
                    densityField[i,j,k] += np.abs(dx*dy*dz)
                    if DEBUG: print("dx,dy,dz: ", dx, dy, dz)
                    if DEBUG: print("dv: ", np.abs(dx*dy*dz))
                    # print("dx, dy, dz: ", dx,dy,dz)
        
    return densityField,x,y,z


""" Plots the Density field for a given axis and value
 
Inputs:
 ------
 densityField : NumPy array
     the calculated density field of the particle distribution

 x: Numpy array
    array of vertices on the x-axis
   
 y: Numpy array
    array of vertices on the y-axis
    
 z: Numpy array
    array of vertices on the z-axis

 axis : string
     The chosen axis
 
 value : NumPy value
     The value at which the slice is taken

"""  

def PlotDensityField2D(densityField, x,y,z, axis, value):
    if axis == 'y':
        extent = [x[0], x[-1], x[0], x[-1]]  # Define the physical coordinates for the plot
        plt.imshow(densityField[:, value, :], origin='lower', cmap='inferno', extent=extent, aspect='auto')
        plt.title(f'Density Slice at Y={y[value]}')
        plt.ylabel('X-axis')
        plt.xlabel('Z-axis')
    elif axis == 'x':
        extent = [x[0], x[-1], x[0], x[-1]]
        plt.imshow(densityField[value, :, :], origin='lower', cmap='inferno', extent=extent, aspect='auto')
        plt.title(f'Density Slice at X={x[value]}')
        plt.ylabel('Y-axis')
        plt.xlabel('Z-axis')
        
    elif axis == 'z':
        extent = [x[0], x[-1], x[0], x[-1]]
        plt.imshow(densityField[:, :, value], origin='lower', cmap='inferno', extent=extent, aspect='auto')
        plt.title(f'Density Slice at Z={z[value]}')
        plt.ylabel('X-axis')
        plt.xlabel('Y-axis')
            
        
    plt.colorbar(label='Density')
    plt.savefig(f"2DDensityPlot{axis}_{value}.png")
    plt.show()

""" Plots the Density field for a given axis and value
 
Inputs:
 ------
 densityField : NumPy array
     the calculated density field of the particle distribution

 x: Numpy array
    array of vertices on the x-axis
   
 y: Numpy array
    array of vertices on the y-axis
    
 z: Numpy array
    array of vertices on the z-axis

 axis : string
     The chosen plane
 
 value1 : NumPy value
     The value at which the slice is taken
     
 value2 : NumPy value
     The value at which the slice is taken
"""  

def PlotDensityField1D(densityField, x,y,z, axis, value1, value2):
    if axis == 'xy':
        plt.plot(z,densityField[value1, value2, :])
        plt.title(f"Density field with x = {value1} and y = {value2}")
        plt.xlabel("Z-axis")
        plt.ylabel("Density")
    elif axis == 'yz':
        plt.plot(x,densityField[:, value1, value2])
        plt.title(f"Density field with y = {value1} and z = {value2}")
        plt.xlabel("X-axis")
        plt.ylabel("Density")
    elif axis == 'zx':
        plt.plot(y,densityField[value2, :, value1])
        plt.title(f"Density field with x = {value2} and z = {value1}")    
        plt.xlabel("Y-axis")
        plt.ylabel("Density")
        
    plt.savefig(f"1DDensityPlot{axis}_{value1}_{value2}.png")
    plt.show()
