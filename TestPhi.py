#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:03:07 2024

@author: vedhasyamuvva
"""
import numpy as np
import matplotlib.pyplot as plt

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
                dist = x[i] **2 + y[j]**2 + z[k]**2
                if (i == 0 and j == 0 and k ==0 ): 
                    phi[i,j,k] = 1
                else: phi[i,j,k] = 1/dist

    return phi
