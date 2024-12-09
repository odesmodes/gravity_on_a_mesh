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

# Parameters
center = [0,0,0]
a = .2
ba = 1
ca = 1
N = 32**3
grid_size = 32
dt = 0.1

particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)
densityField, x,y,z = InitializePoints.CreateDensityField(center, particles, grid_size)

# Assume velocities = v0 are at rest
# Then  v1/2 = dt/2 * F(x) 

# Define an array for velocities that is the same size as the particles array
velocities = np.zeros_like(particles) 

#PLACEHOLDER FUNCTION UNTIL PHI GETS CALCULATED
def test_phi(x):
    return np.sum(x**2, axis = -1)

# Creating a function to calculate the gradient of the potential at position x with a step size for finite difference
# This is because in order to find the velocities we first need to find the acceleration, which we can get from F = ma = -grad(phi)
def gradient(phi, spacing = 1/32):
    grad = np.gradient(phi, spacing)
    return np.array([grad[0],grad[1],grad[2]])

# Creating a function to calculate the force on each particle at position x
# Force is just negative the gradient of the potential

def F(particles, gradient, spacing = 1/32):

    indices = (particles/spacing).astype(int) # Calculate the grid indices for each particle's position and normalize
    indices = np.clip(indices, 0, gradient.shape[:-3][0] - 1) # Constrain the indices so that they don't go outside the boundaries of the gradient field
    forces = -gradient[indices[:, 0], indices[:, 1], indices[:, 2]] # Append gradient values at each axis for the corresponding indices and multiply by a negative 
    return forces

# test
x_test = np.array([[[1,2,3], [1,2,3], [1,2,3]], [[1,2,3], [1,2,3], [1,2,3]], [[1,2,3], [1,2,3], [1,2,3]]])
phi_test = TestPhi.TestPhi2()
grad_test = gradient(phi_test)
forces_test = F(particles, grad_test)
print("Phi", phi_test)
print("Gradient", grad_test)
print("Force", forces_test)

velocities = dt/2 * F(particles)
# THIS IS from the Verlet Method
def xnext(x, v, dt=dt): 
    v_new = v + F(x) * dt
    x_new = x + v_new * dt
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
ax.set_xlim(-15, 15)
ax.set_ylim(-30, 30)
ax.set_zlim(-15, 15)

# Animate
ani = animation.FuncAnimation(fig, update, frames=1, interval=100, blit=False)
plt.show()



#THIS IS THE pseudo code I started with 
#For actual Verlet Method we will start all particles at rest
#STEP 1 how to verlet method
#STEP 2 how to turn it into an animation
# x0 = particles
# v12 = velocities
# dt = 0.1

# x1 = x0 + v12*dt = x0

# #Loop this part:
# vnew = v12 + F(x1)*dt
# xnew = x1 + vnew*dt

# THIS IS THE ACTUAL CODE FROM ABOVE
# for i in range (5):
#     x1, v12 = xnext(x1, v12, dt)
#     print(x1)
