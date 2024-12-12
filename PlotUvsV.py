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
import IntegrationMethods


# Parameters

a = .2
ba = 1
ca = 1
N = 32**3
dt = 0.1
T = 30

center = [0,0,0]
grid_size = 32

particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, 15)
densityField, _,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)

# Assume velocities = v0 are at rest
# Then  v1/2 = dt/2 * F(x) 

# Define an array for velocities that is the same size as the particles array
#print("particles: ", np.shape(particles))
velocities = np.zeros_like(particles) 

phi = PoissonEquation.IsolatedMass(densityField)
velocities = dt/2 * IntegrationMethods.F(particles, IntegrationMethods.gradient(phi))

# THIS IS from the Verlet Method
def xnext(x, v, phi, dt=dt): 
    densityField, _,_,_ = InitializePoints.CreateDensityField(center, x, grid_size)
    phi = PoissonEquation.IsolatedMass(densityField)
    gradArr = IntegrationMethods.gradient(phi)
    

    # VERLET METHOD FROM WIKIPEDIA
    F_old = IntegrationMethods.F(x, gradArr)
    x_new = x + v * dt + 0.5 * F_old * dt**2
    gradArr_new = IntegrationMethods.gradient(phi)
    F_new = IntegrationMethods.F(x_new, gradArr_new)
    v_new = v + 0.5 * (F_old + F_new) * dt

    return x_new, v_new, phi

def U(particles, phi, spacing = 1/32):
    indices = ((particles + 0.5) / spacing).astype(int)
    U = 0
    U += np.linalg.norm(phi[indices])
    return -U

def K(velocities):
    return np.sum(velocities**2)/2
# Initialize and plot initialized particles and create bounds

Uarr = np.empty(T)
Karr = np.empty(T)

for i in range (T):
    Uarr[i] = U(particles, phi)
    Karr[i] = K(velocities)
    particles, velocities, phi = xnext(particles, velocities, phi)

plt.plot(Karr, Uarr,".")
#plt.plot(Karr, -2*Karr)
plt.show()
