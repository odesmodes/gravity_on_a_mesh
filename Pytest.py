#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 20:50:32 2024

@author: vedhasyamuvva
"""
import numpy as np
import InitializePoints
import PoissonEquation
import pytest

#UNIT TESTS

def test_pointAtCenterOfCell():
    """
    Test the behavior of the density field and potential function when a single particle is placed at the center of a grid cell.
    
    This test ensures that:
    1. The density field is non-zero only at the particle's position.
    2. The density field is zero at adjacent cells.
    3. The potential function (phi) is symmetrical around the particle's position.
    4. The potential is consistent in neighboring grid cells.
    """
    center = [0,0,0]
    grid_size = 32
    
    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]
    
    # Place one particle at exactly the middle point of the grid
    particles = np.array([[x[16], y[16], z[16]]])
    
    # Test the density field is nonzero only at that point
    densityField, _, _, _ = InitializePoints.CreateDensityField(center, particles, grid_size)
    assert (densityField[16][16][16] != 0)
    print(densityField[16][16][16])
    assert (densityField[15][15][16] == 0)
    assert (densityField[15][16][15] == 0)
    assert (densityField[16][15][15] == 0)

    # Test that phi is symmetrical around that one point
    phi = PoissonEquation.IsolatedMass(densityField)
    assert (phi[15][16][16] == phi[16][15][16] and 
            phi[16][16][15] == phi[16][15][16])
    
    assert (phi[17][16][16] == phi[16][17][16] and 
            phi[16][16][17] == phi[16][17][16])
    
    assert (phi[15][16][16] == phi[16][17][16])

# Call the test function
test_pointAtCenterOfCell()
