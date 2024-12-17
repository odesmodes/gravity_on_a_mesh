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
    CreateSolidSphereDensityField: Creates theoretical constant density field for a spherical shell
    CreateSphericalShellDensityField: Creates theoretical constant density field for a solid sphere
    PlotDensityField2D: Plots Density field on a 2D colormap
    PlotDensityField1D: Plots Density field along given axes
    PlotTestFields: Plots predetermined test plots for a density field
"""

import numpy as np
import matplotlib.pyplot as plt


def initializeGaussianPoints(center, a, ba, ca, N = 32**3):
    """
    Generates N random points according to a gaussian distribution set to the chosen parameters

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


def plotInitialPoints(arr):
    """
    Plots the points generated with initialize Gaussian Points

    Inputs:
    ------
    arr : NumPy array
        array of N 3D points
    """
    ax = plt.figure().add_subplot(projection='3d')
    arr = arr.T
    ax.scatter(arr[0], arr[1], arr[2], '.', label='initial points')
    
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    plt.title('3D Scatter Plot of Gaussian Distributed Points')
    plt.show()


def CreateDensityField(center, arr, grid_size):
    """
    Calculates the Density field for every point and returns the density field

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
        the density of each cell

    x: Numpy array
        array of vertices on the x-axis

    y: Numpy array
        array of vertices on the y-axis

    z: Numpy array
        array of vertices on the z-axis
    """
    DEBUG = False
    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]

    densityField = np.zeros((grid_size, grid_size, grid_size))

    for ptcl in arr:
        xi = np.searchsorted(x, ptcl[0])
        yi = np.searchsorted(y, ptcl[1])
        zi = np.searchsorted(z, ptcl[2])
        
        if xi == 0: xarr = np.array([xi])
        elif xi == grid_size: xarr = np.array([xi-1])
        else: xarr = np.array([xi-1, xi])
        
        if yi == 0: yarr = np.array([yi])
        elif yi >= grid_size: yarr = np.array([yi-1])
        else: yarr = np.array([yi-1, yi])
        
        if zi == 0: zarr = np.array([zi])
        elif zi >= grid_size: zarr = np.array([zi-1])
        else: zarr = np.array([zi-1, zi])

        L = 1/grid_size

        for i in xarr:
            for j in yarr:
                for k in zarr:
                    point = np.array([x[i],y[j],z[k]])
                    distance = (ptcl-point)
                    
                    dx = max(0, L - np.abs(distance[0]))
                    dy = max(0, L - np.abs(distance[1]))
                    dz = max(0, L - np.abs(distance[2]))
                    
                    densityField[i,j,k] += np.abs(dx*dy*dz)/L**3
    return densityField, x, y, z


def CreateSphericalShellDensityField(grid_size = 32, radius_sq = 0.3**2, tol = 0.02):
    """
    Calculates a theoretical Density field with a spherical shell of constant density

    Inputs:
    ------
    grid_size : NumPy value
        the number of desired points on the mesh grid

    radius_sq : NumPy value
        the squared value of the desired radius of the sphere

    tol : NumPy value
        the tolerance or half of the width of the spherical shell

    Returns:
    -------
    densityField: Numpy array
        the densityField for a spherical shell of constant density
    """
    densityField = np.zeros((grid_size, grid_size, grid_size))
    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True)

    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                dist = x[i]**2 + y[j]**2 + z[k]**2
                if (dist >= radius_sq - tol and dist <= radius_sq + tol):
                    densityField[i,j,k] = 1

    return densityField


def CreateSolidSphereDensityField(grid_size = 32, radius_sq = 0.3**2):
    """
    Calculates a theoretical Density field for a sphere of constant density

    Inputs:
    ------
    grid_size : NumPy value
        the number of desired points on the mesh grid

    radius_sq : NumPy value
        the squared value of the desired radius of the sphere

    Returns:
    -------
    densityField: Numpy array
        the densityField for a solid sphere of constant density
    """
    densityField = np.zeros((grid_size, grid_size, grid_size))
    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True)

    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                dist = x[i]**2 + y[j]**2 + z[k]**2
                if (dist <= radius_sq):
                    densityField[i,j,k] = 1

    return densityField


def PlotDensityField2D(densityField, x, y, z, axis, value):
    """
    Plots the Density field for a given axis and value

    Inputs:
    ------
    densityField : NumPy array
        the calculated density field of the particle distribution

    x : Numpy array
        array of vertices on the x-axis

    y : Numpy array
        array of vertices on the y-axis

    z : Numpy array
        array of vertices on the z-axis

    axis : string
        The chosen axis

    value : NumPy value
        The value at which the slice is taken
    """
    if axis == 'y':
        extent = [x[0], x[-1], x[0], x[-1]]
        plt.imshow(densityField[:,value,:], extent=extent)
        plt.xlabel('x-axis')
        plt.ylabel('z-axis')
    elif axis == 'x':
        extent = [y[0], y[-1], z[0], z[-1]]
        plt.imshow(densityField[value,:,:], extent=extent)
        plt.xlabel('y-axis')
        plt.ylabel('z-axis')
    elif axis == 'z':
        extent = [x[0], x[-1], y[0], y[-1]]
        plt.imshow(densityField[:,:,value], extent=extent)
        plt.xlabel('x-axis')
        plt.ylabel('y-axis')
        
    plt.colorbar(label='Density')
    plt.title(f'Density Field along {axis}-axis at {value}')
    plt.show()

def PlotDensityField1D(densityField, x, y, z, plane, value1, value2):
    """
    Plot a 1D slice of the density field along a specified plane.

    This function generates a 1D plot of the density field by slicing it along the specified plane 
    ('xy', 'yz', or 'zx') at specific grid coordinates defined by `value1` and `value2`.

    Parameters:
    densityField (np.ndarray): A 3D array representing the density field.
    x (np.ndarray): The x-coordinate grid.
    y (np.ndarray): The y-coordinate grid.
    z (np.ndarray): The z-coordinate grid.
    plane (str): The plane to slice along ('xy', 'yz', or 'zx').
    value1 (int): The first index used for slicing in the chosen plane.
    value2 (int): The second index used for slicing in the chosen plane.

    Returns:
    None: Displays the plot of the sliced density field.
    """
    if plane == 'xy':
        plt.plot(z, densityField[value1, value2, :])
        plt.title(f"Density field with x = {round(value1,3)} and y = {round(value2,3)}")
        plt.xlabel("Z-axis")
        plt.ylabel("Density")
    elif plane == 'yz':
        plt.plot(x, densityField[:, value1, value2])
        plt.title(f"Density field with y = {round(value1,3)} and z = {round(value2,3)}")
        plt.xlabel("X-axis")
        plt.ylabel("Density")
    elif plane == 'zx':
        plt.plot(y, densityField[value2, :, value1])
        plt.title(f"Density field with x = {round(value2,3)} and z = {round(value1, 3)}")    
        plt.xlabel("Y-axis")
        plt.ylabel("Density")
        
    #plt.savefig(f"TestDensityPlots/1DDensityPlot{plane}_{round(value1,3)}_{round(value2,3)}.png")
    plt.show()


def PlotTestFields(densityField, x, y, z, grid_size):
    """
    Plot multiple slices and 2D projections of the density field.

    This function generates 1D slices and 2D projections of the density field at the center of the grid
    along the 'xy', 'yz', and 'zx' planes, and then plots 2D projections along the x, y, and z axes.

    Parameters:
    densityField (np.ndarray): A 3D array representing the density field.
    x (np.ndarray): The x-coordinate grid.
    y (np.ndarray): The y-coordinate grid.
    z (np.ndarray): The z-coordinate grid.
    grid_size (int): The size of the grid, used to define the center for plotting.

    Returns:
    None: Displays multiple plots of the density field.
    """
    PlotDensityField1D(densityField, x, y, z, 'xy', int(grid_size/2), int(grid_size/2))
    PlotDensityField1D(densityField, x, y, z, 'yz', int(grid_size/2), int(grid_size/2))
    PlotDensityField1D(densityField, x, y, z, 'zx', int(grid_size/2), int(grid_size/2))
    
    PlotDensityField2D(densityField, x, y, z, "x", int(grid_size/2))
    PlotDensityField2D(densityField, x, y, z, "y", int(grid_size/2))
    PlotDensityField2D(densityField, x, y, z, "z", int(grid_size/2))
