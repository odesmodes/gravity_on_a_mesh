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
    
""" Calculates the Density field for a single point in the array
 
Inputs:
 ------
 point : NumPy array
     The point at which the density will be calculated

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
 dens: Numpy value
    the density of of the cell around "point"

"""    
def calculateDensityOfPoint(point, a, ba, ca, arr):
    #Theoretically these values should be divided by 32
    La = a
    Lb = ba*a
    Lc = ca*a
    #print(La,Lb,Lc)
    distance = np.abs(arr-point)
    #print(arr)
    #print(distance)
    dens = 0
    
    for i in range (np.shape(distance)[0]):
        if (distance[i][0]<= (La/2) and distance[i][1] <= (Lb/2) and distance[i][2] <= (Lc/2)):
            dens += ((La - distance[i][0])*(Lb - distance[i][1])*(Lc - distance[i][2]))
    return dens

""" Calculates the Density field for every point and plot some results
 
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
def PlotDensityField(center, a, ba, ca, arr):
    #Create a a NxNxN Grid of points to calculate the density for 
    x = np.arange(32) * a + center[0]
    y = np.arange(32) * ba * a + center[1]
    z = np.arange(32) * ca * a + center[2]

    # Create the meshgrid
    X, Y, Z = np.meshgrid(x, y, z)

    # Stack X, Y, Z along a new axis to create a grid of grids
    #points = np.stack([X, Y, Z], axis=-1)
    grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
    
    #Get values for each point
    densityField = [calculateDensityOfPoint(point, a, ba, ca, arr) for point in grid_points]
    
    #Plot along each of the axes
    plt.plot(grid_points[:,16,16],densityField[:,16,16], '.')
    plt.show()
    
    
    return densityField
