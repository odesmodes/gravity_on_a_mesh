#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:42:15 2024

@author: vedhasyamuvva
"""

# Section 2 Part 1: Distributing particles according to a multi-variate Gaussian with a chosen center, semimajor axis a, and axis ratios b/a and c/a. 

import numpy as np
import matplotlib.pyplot as plt


def initializeGaussianPoints(center, a, ba, ca, N = 32**3):
    mean = np.array(center)
    cov = np.diag(np.array([a**2, (a*ba)**2, (a*ca)**2]))
    
    return np.random.multivariate_normal(mean, cov, N)

def plotInitialPoints(arr):
    
    ax = plt.figure().add_subplot(projection='3d')
    arr = arr.T
    ax.scatter(arr[0], arr[1], arr[2], '.', label='initial points')
    
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    plt.title('3D Scatter Plot of Gaussian Distributed Points')
    
    
    plt.show()
