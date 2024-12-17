#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:03:07 2024

@author: vedhasyamuvva
"""
import numpy as np
import matplotlib.pyplot as plt
import InitializePoints

def TestPhi1(grid_size = 32, radius_sq = 0.3**2):
    """
    Generate a potential field for a spherical region with non-zero values inside a defined radius.

    This function calculates a potential field where the potential is equal to the square of 
    the distance from the origin for points inside a defined radius and zero outside.

    Parameters:
    grid_size (int): The size of the grid for the potential field.
    radius_sq (float): The square of the radius inside which the potential is calculated.

    Returns:
    np.ndarray: A 3D array representing the potential field.
    """
    phi = np.zeros((grid_size, grid_size, grid_size))
    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    
    for i in range (grid_size):
        for j in range (grid_size):
            for k in range (grid_size):
                dist = x[i] **2 + y[j]**2 + z[k]**2
                if (dist <= radius_sq):
                    phi[i,j,k] = dist

    return phi


def TestPhi2(grid_size = 32):
    """
    Generate a potential field for a point source located at the center of the grid.

    This function computes the inverse square potential of a point mass located at the center of the grid,
    with a specific value at the origin to avoid division by zero.

    Parameters:
    grid_size (int): The size of the grid for the potential field.

    Returns:
    np.ndarray: A 3D array representing the potential field with the inverse-square form.
    """
    phi = np.zeros((grid_size, grid_size, grid_size))
    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
    
    for i in range (grid_size):
        for j in range (grid_size):
            for k in range (grid_size):
                dist_sq = float(x[i] **2 + y[j]**2 + z[k]**2)
                if (dist_sq == 0): 
                    phi[i,j,k] = -400
                else: 
                    phi[i,j,k] = -1/dist_sq

    return phi 


def TestPhi3(grid_size, mass=1.0):
    """
    Generate a 3D gravitational potential field for a single mass centered on a grid.

    This function computes the gravitational potential at each grid point due to a single mass 
    located at the center of the grid, using the formula for gravitational potential.

    Parameters:
    grid_size (int): The size of the grid for the potential field.
    mass (float): The mass at the center of the grid.

    Returns:
    np.ndarray: A 3D array representing the gravitational potential field.
    """
    # Determine the center of the grid
    center = np.array([grid_size // 2] * 3)

    # Create a 3D grid of coordinates
    x, y, z = np.meshgrid(
        np.arange(grid_size),
        np.arange(grid_size),
        np.arange(grid_size),
        indexing='ij'
    )

    # Calculate the distance from the center for all grid points
    r_squared = (x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2

    # Avoid division by zero at the center
    r_squared[r_squared == 0] = np.finfo(float).eps

    # Compute the gravitational potential
    phi = -mass / np.sqrt(r_squared)

    return phi


def TestPhi4(grid_size=32, mass=100, softening=1e-2):    
    """
    Generate a 3D potential field for a single mass at the center with a softening factor to avoid singularities.

    This function computes the potential for a point mass at the center of the grid, 
    with a softening parameter to prevent division by zero at the center.

    Parameters:
    grid_size (int): The size of the grid for the potential field.
    mass (float): The mass of the point source located at the center of the grid.
    softening (float): A small value added to the distance to prevent singularities at the center.

    Returns:
    np.ndarray: A 3D array representing the potential field with softening.
    """
    # Create a grid of coordinates
    x = np.linspace(-0.5, 0.5, grid_size)
    y = np.linspace(-0.5, 0.5, grid_size)
    z = np.linspace(-0.5, 0.5, grid_size)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    # Calculate the distance to the center of the grid
    center = (0.0, 0.0, 0.0)  # Center of the grid in normalized space
    distance = np.sqrt((X - center[0])**2 + (Y - center[1])**2 + (Z - center[2])**2 + softening**2)

    # Compute the potential field
    potential = -mass / distance

    return potential


def gradient(phi, spacing = 1/32):
    """
    Calculate the gradient of the potential field.

    This function calculates the gradient of the potential field in three dimensions using
    numerical differentiation, assuming a uniform grid spacing.

    Parameters:
    phi (np.ndarray): The 3D potential field for which the gradient is to be computed.
    spacing (float): The spacing between grid points in each dimension.

    Returns:
    np.ndarray: A 3D array representing the gradient of the potential field.
    """
    grad = np.gradient(phi, spacing)
    return np.array([grad[0], grad[1], grad[2]])


"""
grid_size = 32

x = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
y = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
z = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
phi = TestPhi2()
InitializePoints.PlotTestFields(phi, x,y,z,grid_size)
"""
