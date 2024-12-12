#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:38:05 2024

@author: vedhasyamuvva
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import InitializePoints
import TestPhi
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
def xnext(x, v, dt=dt): 
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
    return x_new, v_new

# Update for each frame, particles and velocities get updated
def update(frame):
    global particles, velocities
    particles, velocities = xnext(particles, velocities, dt)
    scat._offsets3d = (particles[:, 0], particles[:, 1], particles[:, 2])
    return scat,

# Initialize and plot initialized particles and create bounds
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
scat = ax.scatter(particles[:, 0], particles[:, 1], particles[:, 2], label='Particles')
ax.set_xlim(-0.6, 0.6)
ax.set_ylim(-0.6, 0.6)
ax.set_zlim(-0.6, 0.6)

ax.set_title("Exploring The Physics Of Gravitational Collapse")
# Animate
ani = animation.FuncAnimation(fig, update, frames=5, interval=100, blit=False)
ani.save('grav_collapse.mp4', fps=20, extra_args=['-vcodec', 'libx264'])
plt.show()


