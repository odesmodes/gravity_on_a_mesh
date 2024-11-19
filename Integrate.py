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

# Parameters
center = [0,0,0]
a = 4
ba = 2
ca = 1
N = 32**3
grid_size = 32


particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)
# Assume velocities = v_1/2 are at rest
velocities = np.zeros_like(particles)
# Since particles start at rest, x0 = x1

dt = 0.1

# THIS NEEDS TO BE UPDATED ONCE WE HAVE FINISHED PROBLEM 3
def F(x):
    # Negative gradient of a potential function ONCE PROBLEM 3 is done
    return -2 * x + 3

# THIS IS from the Verlet Method
def xnext(x, v, dt=0.1): 
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
ani = animation.FuncAnimation(fig, update, frames=5, interval=1, blit=False)
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