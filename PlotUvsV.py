#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:22:15 2024

@author: vedhasyamuvva
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:38:05 2024

@author: vedhasyamuvva
"""
import matplotlib.pyplot as plt
import numpy as np
import InitializePoints
import PoissonEquation


# Creating a function to calculate the gradient of the potential at position x with a step size for finite difference
# This is because in order to find the velocities we first need to find the acceleration, which we can get from F = ma = -grad(phi)
def gradient(phi, spacing = 1/32):
    grad = np.gradient(phi, spacing)
    return np.array([grad[0],grad[1],grad[2]])

# Creating a function to calculate the force on each particle at position x
# Force is just negative the gradient of the potential

def F(particles, gradArr, spacing=1/32):
    # Calculate the grid indices for each particle's position
    indices = ((particles + 0.5) / spacing).astype(int)
    
    # Initialize forces array
    forces = np.zeros_like(particles)
    
    # Check if indices are within valid bounds
    in_bounds = (indices[:, 0] >= 0) & (indices[:, 0] < gradArr.shape[1]) & \
                (indices[:, 1] >= 0) & (indices[:, 1] < gradArr.shape[2]) & \
                (indices[:, 2] >= 0) & (indices[:, 2] < gradArr.shape[3])
    
    # Apply forces only to particles within bounds
    valid_indices = indices[in_bounds]
    forces[in_bounds] = gradArr[:, valid_indices[:, 0], valid_indices[:, 1], valid_indices[:, 2]].T
    
    return -forces

# Parameters
center = [0,0,0]
a = .2
ba = 1
ca = 1
N = 32**3
grid_size = 32
dt = 10

particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, 15)
densityField, _,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)

# Assume velocities = v0 are at rest
# Then  v1/2 = dt/2 * F(x) 

# Define an array for velocities that is the same size as the particles array
#print("particles: ", np.shape(particles))
velocities = np.zeros_like(particles) 

phi = PoissonEquation.IsolatedMass(densityField)
velocities = dt/2 * F(particles, gradient(phi))

# THIS IS from the Verlet Method
def xnext(x, v, phi, dt=dt): 
    #RECALCULATE PHI
    densityField, _,_,_ = InitializePoints.CreateDensityField(center, x, grid_size)
    phi = PoissonEquation.IsolatedMass(densityField)
    gradArr = gradient(phi)
    
    # Calculate gradArr with new distribution

    # VERLET METHOD FROM WIKIPEDIA
    F_old = F(x, gradArr)
    x_new = x + v * dt + 0.5 * F_old * dt**2
    gradArr_new = gradient(phi)
    F_new = F(x_new, gradArr_new)
    v_new = v + 0.5 * (F_old + F_new) * dt

    #v_new = v + F(x, gradArr) * dt/2
    #print("x: ", np.shape(x), np.shape(v_new))
    #x_new = x + v_new * dt
    return x_new, v_new, phi

def U(particles, phi, spacing = 1/32):
    indices = ((particles + 0.5) / spacing).astype(int)
    U = 0
    U += np.linalg.norm(phi[indices])
    return U

def K(velocities):
    return np.sum(velocities**2)
# Initialize and plot initialized particles and create bounds

T = 10
Uarr = np.empty(T)
Karr = np.empty(T)

for i in range (T):
    Uarr[i] = U(particles, phi)
    Karr[i] = K(velocities)
    particles, velocities, phi = xnext(particles, velocities, phi)

plt.plot(Karr, Uarr,".")
plt.plot(Karr, 2*Karr+0.007)
plt.show()
