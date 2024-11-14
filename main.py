#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:40:32 2024

@author: vedhasyamuvva
"""

import InitializePoints
import numpy as np

center = [0,0,0]
a = 1
ba = 2
ca = 3
N = 32**3


arr = InitializePoints.initializeGaussianPoints(center, a, ba, ca)

InitializePoints.plotInitialPoints(arr)

testpoint = np.array([0,0,0])
#print(InitializePoints.calculateDensityOfPoint(testpoint,a,ba,ca, arr))

InitializePoints.PlotDensityField(center, a, ba, ca, arr, 'x', 32)
