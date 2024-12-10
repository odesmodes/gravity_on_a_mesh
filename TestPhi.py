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
                if ( dist <= radius_sq):
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
