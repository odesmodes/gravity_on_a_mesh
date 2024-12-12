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


velocities = np.zeros_like(particles) 

phi = PoissonEquation.IsolatedMass(densityField)
velocities = dt/2 * IntegrationMethods.F(particles, IntegrationMethods.gradient(phi))


# Update for each frame, particles and velocities get updated
def update(frame):
    global particles, velocities
    particles, velocities = IntegrationMethods.xnext(particles, velocities, dt)
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
#ani.save('grav_collapse.mp4', fps=20, extra_args=['-vcodec', 'libx264'])
plt.show()


