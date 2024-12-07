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
velocities = np.zeros_like(particles)

#PLACEHOLDER FUNCTION UNTIL PHI GETS CALCULATED
def phi(x):
    return np.sum(x**2)

# Creating a function to calculate the gradient of the potential at position x with a step size for finite difference

def gradient(phi, x, h=1e-5):

	dim = len(x)
	grad = np.zeros(dim)
	
	for i in range(dim): 
		x_f = x.copy() # makes copies of x to avoid messing up our og x array
		x_i = x.copy()
		x_f[i] += h # add increments of h i.e. dx
		x_i[i] -= h
		grad[i] = (phi(x_f)-phi(x_i))/(2*h) # def of gradient
	return grad

# Creating a function to calculate the force on each particle at position x
# Force is just negative the gradient of the potential

def F(x):
	forces = np.zeros_like(x)
	
	for i, x0  in enumerate(x):
		forces[i] = -gradient(phi,x0)
	return forces

# test
x_test = np.array([1.0,2.0,3.0])
grad_test = gradient(phi, x_test)
forces_test = F(np.array([x_test]))
print("Gradient at", x_test, "is", grad_test)
print("Force at", x_test, "is", forces_test)


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
