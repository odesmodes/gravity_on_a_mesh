#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:40:32 2024

@author: vedhasyamuvva
"""

import InitializePoints
import numpy as np

center = [0,0,0]
a = 4
ba = 1
ca = 1
N = 32**3
grid_size = 32



arr = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)
# testpoint = np.array( [0, 0, 0])
# testpoint2 = np.array([0, 0,  .5])
# testpoint3 = np.array([0, 0,  -1.5])

# arr = np.array([testpoint, testpoint2, testpoint3])

InitializePoints.plotInitialPoints(arr)

densityField, x,y,z = InitializePoints.CreateDensityField(center, a, ba, ca, arr, grid_size)

#print("FinalDensityField: ", densityField)    
            
InitializePoints.PlotDensityField1D(densityField, x,y,z,'xy', int(grid_size/2),int(grid_size/2))
InitializePoints.PlotDensityField1D(densityField, x,y,z,'yz', int(grid_size/2),int(grid_size/2))
InitializePoints.PlotDensityField1D(densityField, x,y,z,'zx', int(grid_size/2),int(grid_size/2))

InitializePoints.PlotDensityField2D(densityField, x,y,z, "x", int(grid_size/2))
InitializePoints.PlotDensityField2D(densityField, x,y,z, "y", int(grid_size/2))
InitializePoints.PlotDensityField2D(densityField, x,y,z, "z", int(grid_size/2))
