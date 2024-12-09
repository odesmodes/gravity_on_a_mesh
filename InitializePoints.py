#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:42:15 2024

@author: vedhasyamuvva
"""

"""
______________________________________________________________________________
Includes the methods to create initial random points and calculate the Density field
______________________________________________________________________________

methods:
    initializeGaussianPoints: Creates a Gaussian Distribution of points
    plotInitialPoints: plots particles on a 3D scatter plot
    CreateDensityField: Calculates density fields for given particle distribution
    PlotDensityField2D: Plots Density field on a 2D colormap
    PlotDensityField1D: Plots Density field along given axes
    PlotTestFields: Plots predetermined test plots for a density field
"""

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
    points = np.empty((N,3))
    hasSpace = True
    i = 0
    while hasSpace:
        point = np.random.multivariate_normal(mean, cov)
        if (point[0] >= -0.5 and point[0] <= 0.5 and 
            point[1] >= -0.5 and point[1] <= 0.5 and 
            point[2] >= -0.5 and point[2] <= 0.5): 
            points[i] = point
            i += 1
            if (i >= N): hasSpace = False
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
def CreateDensityField(center, arr, grid_size):
    #FOR DEBUGGING PURPOSES:
    DEBUG = False

    #Create a NxNxN mesh grid to be overlayed over particles
    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]

    
    #Initialize Density Field
    densityField = np.zeros((grid_size, grid_size, grid_size))
    

    # Calculate the addition to the density field from each individual particle
    for ptcl in arr:
        
        #Find the grid point closest to the particle
        
        # print("particle: ", ptcl)
        xi = np.searchsorted(x, ptcl[0])
        yi = np.searchsorted(y, ptcl[1])
        zi = np.searchsorted(z, ptcl[2])
        if DEBUG: print("xi,yi,zi: ", xi,yi,zi)
        
        #Check if the particle is out of bounds
        if xi == 0: xarr = np.array([xi])
        elif xi == grid_size: xarr = np.array([xi-1])
        else: xarr = np.array([xi-1,xi])
        if DEBUG: print("xarr: ", xarr)
        
        if yi == 0: yarr = np.array([yi])
        elif yi >= grid_size: yarr = np.array([yi-1])
        else: yarr = np.array([yi-1,yi])
        if DEBUG: print("yarr: ", yarr)
        
        if zi == 0: zarr = np.array([zi])
        elif zi >= grid_size: zarr = np.array([zi-1])
        else: zarr = np.array([zi-1,zi])
        if DEBUG: print("zarr: ", zarr)
        
        # L is the size of each cell in the mesh grid
        L = 1/grid_size
        
        #Calculate the density addition at each vertex determined above
        for i in xarr:
            for j in yarr:
                for k in zarr:
                    if DEBUG: print("ijk: ", i, j, k)
                    
                    point = np.array([x[i],y[j],z[k]])
                    
                    if DEBUG: print("ptcl: ", ptcl)
                    if DEBUG: print("point: ", point)
                    
                    distance = (ptcl-point)
                    
                    if DEBUG: print("distance: ", distance)
                    
                    #Overlap in cubical cloud in each direction
                    dx = max(0, L - np.abs(distance[0]))
                    dy = max(0, L - np.abs(distance[1]))
                    dz = max(0, L - np.abs(distance[2]))
                    
                    #The addition to the density at the vertex from the particle
                    densityField[i,j,k] += np.abs(dx*dy*dz)
                    
                    if DEBUG: print("dx,dy,dz: ", dx, dy, dz)
                    if DEBUG: print("dv: ", np.abs(dx*dy*dz))
        
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
        plt.title(f'Density Slice at Y={round(y[value],3)}')
        plt.ylabel('X-axis')
        plt.xlabel('Z-axis')
    elif axis == 'x':
        extent = [x[0], x[-1], x[0], x[-1]]
        plt.imshow(densityField[value, :, :], origin='lower', cmap='inferno', extent=extent, aspect='auto')
        plt.title(f'Density Slice at X={round(x[value],3)}')
        plt.ylabel('Y-axis')
        plt.xlabel('Z-axis')
        
    elif axis == 'z':
        extent = [x[0], x[-1], x[0], x[-1]]
        plt.imshow(densityField[:, :, value], origin='lower', cmap='inferno', extent=extent, aspect='auto')
        plt.title(f'Density Slice at Z={round(z[value],3)}')
        plt.ylabel('X-axis')
        plt.xlabel('Y-axis')
            
        
    plt.colorbar(label='Density')
    #plt.savefig(f"TestDensityPlots/2DDensityPlot{axis}_{round(value,3)}.png")
    plt.show()




""" Plots the Density field for a given plane and value
 
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

 plane : string
     The chosen plane
 
 value1 : NumPy value
     The value at which the slice is taken
     
 value2 : NumPy value
     The value at which the slice is taken
"""  

def PlotDensityField1D(densityField, x,y,z, plane, value1, value2):
    if plane == 'xy':
        plt.plot(z,densityField[value1, value2, :])
        plt.title(f"Density field with x = {round(value1,3)} and y = {round(value2,3)}")
        plt.xlabel("Z-axis")
        plt.ylabel("Density")
    elif plane == 'yz':
        plt.plot(x,densityField[:, value1, value2])
        plt.title(f"Density field with y = {round(value1,3)} and z = {round(value2,3)}")
        plt.xlabel("X-axis")
        plt.ylabel("Density")
    elif plane == 'zx':
        plt.plot(y,densityField[value2, :, value1])
        plt.title(f"Density field with x = {round(value2,3)} and z = {round(value1, 3)}")    
        plt.xlabel("Y-axis")
        plt.ylabel("Density")
        
    #plt.savefig(f"TestDensityPlots/1DDensityPlot{plane}_{round(value1,3)}_{round(value2,3)}.png")
    plt.show()


""" Creates 1D and 2D plots of the density fields around the center point
 
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

 grid_size : NumPy value
     The size of the grid
"""  

def PlotTestFields(densityField, x,y,z, grid_size):
    PlotDensityField1D(densityField, x,y,z,'xy', int(grid_size/2),int(grid_size/2))
    PlotDensityField1D(densityField, x,y,z,'yz', int(grid_size/2),int(grid_size/2))
    PlotDensityField1D(densityField, x,y,z,'zx', int(grid_size/2),int(grid_size/2))
    
    PlotDensityField2D(densityField, x,y,z, "x", int(grid_size/2))
    PlotDensityField2D(densityField, x,y,z, "y", int(grid_size/2))
    PlotDensityField2D(densityField, x,y,z, "z", int(grid_size/2))