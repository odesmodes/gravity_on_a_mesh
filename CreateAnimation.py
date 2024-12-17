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
import IntegrationMethods


# Parameters for particle distribution and simulation
center = [0,0,0]  # Center of the particle distribution
a = .2             # Semi-major axis for Gaussian distribution
ba = 1             # Axis ratio for the second axis
ca = 1             # Axis ratio for the third axis
usedN = 15             # Number of particles to be generated
grid_size = 32     # Grid size for density field
dt = .01           # Time step for the simulation

# Initialization of particles and density field
numStep = 1
particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)
densityField, _,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)

# Initialize particle velocities
velocities = np.zeros_like(particles) 
# velocities = IntegrationMethods.NetAngularMomentum(particles)

# Calculate the gravitational potential (phi) and update velocities
phi = PoissonEquation.IsolatedMass(densityField)
velocities = velocities + dt/2 * IntegrationMethods.F(particles, IntegrationMethods.gradient(phi))

# Update function for each animation frame, where particle positions and velocities are updated
def update(frame):
    """
    Updates the particle positions and velocities for each frame in the animation.

    Parameters:
    -----------
    frame : int
        The current frame number in the animation sequence.

    Returns:
    --------
    tuple
        A tuple containing the scatter object which updates the plot in each frame.
    """
    global particles, velocities, numStep, dt
    particles, velocities = IntegrationMethods.xnext(particles, velocities, dt)
    scat._offsets3d = (particles[:, 0], particles[:, 1], particles[:, 2])
    numStep += 1
    if (numStep % 10 == 0): dt = dt * 0.9  # Decay time step every 10 steps
    return scat,

# Initialize the figure, plot particles, and set plot bounds
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
scat = ax.scatter(particles[:, 0], particles[:, 1], particles[:, 2], label='Particles')
ax.set_xlim(-0.6, 0.6)
ax.set_ylim(-0.6, 0.6)
ax.set_zlim(-0.6, 0.6)

ax.set_title("Exploring The Physics Of Gravitational Collapse")

# Animate the particle system over 50 frames with an interval of 100ms
ani = animation.FuncAnimation(fig, update, frames=50, interval=100, blit=False)

# Display the animation
# ani.save('Vgrav_collapse.mp4', fps=5)
plt.show()
