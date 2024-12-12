#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 20:50:32 2024

@author: vedhasyamuvva
"""
import numpy as np
import InitializePoints
import pytest

#UNIT TESTS
def test_pointAtCenterOfCell():
    center = [0,0,0]
    grid_size = 32
    
    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]
    
    particles = np.array([[x[16],y[16],z[16]]])
    
    densityField,_,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)
    assert (densityField[16][16][16] == 1)
    print(densityField[16][16][16])
    assert (densityField[15][15][16] == 0)
    assert (densityField[15][16][15] == 0)
    assert (densityField[16][15][15] == 0)

test_pointAtCenterOfCell()
