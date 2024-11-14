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

 Returns:
 -------
 densityField: Numpy array
    the density of of the each cell arranged in a somewhat convoluted fashion that needs to be fixed

"""  
def CreateDensityField(center, a, ba, ca, arr):
    N = 32
    #Create a a NxNxN Grid of points to calculate the density for 
    x = np.linspace(center[0] - a, center[0] + a, N)
    y = np.linspace(center[1] - ba * a, center[1] + ba * a, N)
    z = np.linspace(center[2] - ca * a, center[2] + ca * a, N)

    #Initialize Density Field
    densityField = np.zeros((N, N, N))

    La = a/N
    Lb = ba*a/N
    Lc = ca*a/N
    
    # Calculate the addition to the density field from each individual particle
    for ptcl in arr:
        #Find the grid point closest to the particle
        xi = np.searchsorted(x, ptcl[0])
        yi = np.searchsorted(y, ptcl[1])
        zi = np.searchsorted(z, ptcl[2])
        
        
        #Check to see which vertices of the box will matter
        if xi == 0: xarr = np.array([xi])
        elif xi == N: xarr = np.array([xi-1])
        else: xarr = np.array([xi-1,xi])
        
        if yi == 0: yarr = np.array([yi])
        elif yi >= N: yarr = np.array([yi-1])
        else: yarr = np.array([yi-1,yi])
        
        if zi == 0: zarr = np.array([zi])
        elif zi >= N: zarr = np.array([zi-1])
        else: zarr = np.array([zi-1,zi])
        
        #Check the vertices of the box around it
        for i in xarr:
            for j in yarr:
                for k in zarr:
                    point = np.array([x[i],y[i],z[i]])
                    distance = (ptcl-point)
                    dx = max(0, La - distance[0])
                    dy = max(0, La - distance[0])
                    dz = max(0, La - distance[0])
                    
                    densityField[i,j,k] += dx*dy*dz
        
    return densityField


""" Plots the Density field for a given axis and value
 
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

 axis : string
     The chosen axis
 
 value : NumPy value
     The value at which the slice is taken

"""  

def PlotDensityField(center, a, ba, ca, arr, axis, value):
    densityField = CreateDensityField(center, a, ba, ca, arr)
    
    if axis == 'z':
        plt.imshow(densityField[:, :, value], origin='lower', cmap='hot')
        plt.title(f'Density Slice at Z={value}')
    elif axis == 'y':
        plt.imshow(densityField[:, value, :], origin='lower', cmap='hot')
        plt.title(f'Density Slice at Y={value}')
    elif axis == 'x':
        plt.imshow(densityField[value, :, :], origin='lower', cmap='hot')
        plt.title(f'Density Slice at X={value}')
    plt.colorbar(label='Density')
    
    plt.show()
