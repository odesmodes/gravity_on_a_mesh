#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:40:32 2024

@author: vedhasyamuvva
"""

import InitializePoints
import numpy as np

TEST1 = False #Toggle for each test
TEST2 = False #Toggle for each test
TEST3 = False #Toggle for each test
TEST4 = False  #Toggle for each test
TEST5 = True  #Toggle for each test
TEST6 = False  #Toggle for each test

center = [0,0,0]
grid_size = 32


if TEST1:
    particles = np.array([[0,0,0]])


if TEST2:
    particles = np.array([[0,0,0],[0.3,0,0], [0,0.3,0], [0,0,0.3]])


if TEST3:
    a = .2
    ba = 1
    ca = 1
    N = 32**3
    
    particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)
    
if TEST4:
    a = .15
    ba = .5
    ca = 2
    N = 32**3
    
    particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)

if TEST5:
    a = .1
    ba = 1
    ca = .3
    N = 32**3
    
    particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)

if TEST6:
    a = .1
    ba = 2
    ca = .6
    N = 32**3
    
    particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)


InitializePoints.plotInitialPoints(particles)

densityField, x,y,z = InitializePoints.CreateDensityField(center, particles, grid_size)

InitializePoints.PlotTestFields(densityField, x, y, z, grid_size)

