#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:40:32 2024

@author: vedhasyamuvva
"""

import InitializePoints
import numpy as np

# Toggle variables to switch between different tests
TEST1 = True  # Toggle for Test 1
TEST2 = False  # Toggle for Test 2
TEST3 = False  # Toggle for Test 3
TEST4 = False   # Toggle for Test 4
TEST5 = False  # Toggle for Test 5
TEST6 = False  # Toggle for Test 6

# Center coordinates for point distribution
center = [0,0,0]
# Grid size for density field
grid_size = 32

if TEST1:
    """
    Test 1: Creates a single particle at the origin.
    """
    particles = np.array([[0,0,0]])

if TEST2:
    """
    Test 2: Creates four particles at various positions along the axes.
    """
    particles = np.array([[0,0,0],[0.3,0,0], [0,0.3,0], [0,0,0.3]])

if TEST3:
    """
    Test 3: Generates particles with a Gaussian distribution using the specified parameters.
    """
    a = .2
    ba = 1
    ca = 1
    N = 32**3
    particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)

if TEST4:
    """
    Test 4: Generates particles with a Gaussian distribution with different parameters.
    """
    a = .15
    ba = .5
    ca = 2
    N = 32**3
    particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)

if TEST5:
    """
    Test 5: Generates particles with a Gaussian distribution using another set of parameters.
    """
    a = .1
    ba = 1
    ca = .3
    N = 32**3
    particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)

if TEST6:
    """
    Test 6: Generates particles with a Gaussian distribution using another set of parameters.
    """
    a = .1
    ba = 2
    ca = .6
    N = 32**3
    particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)

# Plot the initial particles
InitializePoints.plotInitialPoints(particles)

# Calculate the density field from the generated particles
densityField, x, y, z = InitializePoints.CreateDensityField(center, particles, grid_size)

# Plot test fields to visualize the density distribution
InitializePoints.PlotTestFields(densityField, x, y, z, grid_size)
