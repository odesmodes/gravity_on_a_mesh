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
                else: phi[i,j,k] = -1/dist_sq

    return phi 

def TestPhi3(grid_size, mass=1.0):
    """
    Generate a 3D gravitational potential field for a single mass centered on a grid.

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

# phi = TestPhi3(32)
# print(phi)

import numpy as np

def TestPhi4(grid_size=32, mass=100, softening=1e-3):    

    # Also single mass in center potential take two 
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

grid_size = 32
potential = TestPhi4()

# Check potential at the center
print("Potential at the center:", potential[grid_size // 2, grid_size // 2, grid_size // 2])

# Visualize a slice of the potential field 
#plt.imshow(potential[:, :, grid_size // 2], extent=(-0.5, 0.5, -0.5, 0.5), origin='lower')
#plt.colorbar(label="Potential")
#plt.title("Potential Slice at Z = 0")
#plt.xlabel("X")
#plt.ylabel("Y")
#plt.show()


def gradient(phi, spacing = 1/32):
    grad = np.gradient(phi, spacing)
    return np.array([grad[0],grad[1],grad[2]])

"""
grid_size = 32

x = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
y = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
z = np.linspace(-0.5, 0.5, grid_size, endpoint=True)
phi = TestPhi2()
InitializePoints.PlotTestFields(phi, x,y,z,grid_size)
"""
