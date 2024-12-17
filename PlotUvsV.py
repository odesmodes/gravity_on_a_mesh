#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:22:15 2024

@author: vedhasyamuvva
"""

import matplotlib.pyplot as plt
import numpy as np
import InitializePoints
import PoissonEquation
import IntegrationMethods
from scipy.interpolate import RegularGridInterpolator

# Parameters

a = .2
ba = 1
ca = 1
N = 32**3
dt = 0.00001
T = 10000

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
    """
    Calculates the next position and velocity of particles using the Verlet integration method.

    Inputs:
    ------
    x : NumPy array
        Current positions of the particles.

    v : NumPy array
        Current velocities of the particles.

    phi : NumPy array
        Potential field to calculate the gradient.

    dt : float, optional
        Time step for the integration (default is dt).

    Returns:
    -------
    x_new : NumPy array
        Updated positions of the particles.

    v_new : NumPy array
        Updated velocities of the particles.

    phi : NumPy array
        Updated potential field after recalculating.
    """
    densityField, _, _, _ = InitializePoints.CreateDensityField(center, x, grid_size)
    phi = PoissonEquation.IsolatedMass(densityField)
    gradArr = IntegrationMethods.gradient(phi)
    
    # Verlet method from Wikipedia
    F_old = IntegrationMethods.F(x, gradArr)
    x_new = x + v * dt + 0.5 * F_old * dt**2
    gradArr_new = IntegrationMethods.gradient(phi)
    F_new = IntegrationMethods.F(x_new, gradArr_new)
    v_new = v + 0.5 * (F_old + F_new) * dt

    return x_new, v_new, phi

def U(particles, phi, grid_spacing=1/32, grid_origin=-0.5):
    """
    Calculates the potential energy of the system based on the particle positions and potential field.

    Inputs:
    ------
    particles : NumPy array
        The positions of the particles.

    phi : NumPy array
        The potential field at each point in the grid.

    grid_spacing : float, optional
        The spacing between grid points (default is 1/32).

    grid_origin : float, optional
        The origin of the grid (default is -0.5).

    Returns:
    -------
    potential_energy : float
        The total potential energy of the system.
    """
    # Create coordinates for the grid
    grid_points = np.linspace(grid_origin, grid_origin + grid_spacing * phi.shape[0], phi.shape[0])
    interp_func = RegularGridInterpolator((grid_points, grid_points, grid_points), phi, bounds_error=False, fill_value=0)
    
    # Interpolate phi at particle positions
    potential_values = interp_func(particles)
    
    # Sum potential energy contributions
    return np.sum(potential_values)

def K(velocities):
    """
    Calculates the kinetic energy of the system based on the velocities of the particles.

    Inputs:
    ------
    velocities : NumPy array
        The velocities of the particles.

    Returns:
    -------
    kinetic_energy : float
        The total kinetic energy of the system.
    """
    return np.sum(velocities**2) / 2

# Initialize and plot initialized particles and create bounds

Uarr = np.empty(T)
Karr = np.empty(T)

for i in range(T):
    Uarr[i] = U(particles, phi)
    Karr[i] = K(velocities)
    particles, velocities, phi = xnext(particles, velocities, phi)

plt.plot(Karr, Uarr, ".")
a, b = np.polyfit(Karr, Uarr, 1)

plt.plot(Karr, a * Karr + b, label=f"Potential = {a.round(3)} Kinetic + {b.round(2)}")
plt.xlabel("Kinetic Energy (J)")
plt.ylabel("Potential Energy (J)")
plt.legend()
plt.savefig("GenericUvsV.png")
plt.show()

# CODE TO TEST CONSERVATION OF ENERGY
arr = Uarr + Karr
tarr = np.arange(0, T * dt, dt)
plt.plot(tarr, arr, ".", label="total energy")
plt.xlabel("time")
plt.ylabel("total energy")
plt.legend()
plt.show()
print(dt)
