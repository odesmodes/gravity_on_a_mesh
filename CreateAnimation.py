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


"""
______________________________________________________________________________
Run this file to create an animation Exploring the Gravitational Collapse
______________________________________________________________________________

Adjustable Paramters:
    a: Numpy Value
        semi major axis
        
    ba: Numpy Value 
        first semi minor to semi major axis ratio
    
    ca: Numpy Value
        second semi minor to semi major axis ratio
        
    N: Number Value
        number of particles in the simulation
    
    dt: Numpy Value
        time step
    
"""

center = [0,0,0] #Define center of the diagram (usually set to [0,0,0])
a = .2   #Value of semimajor axis
ba = 1   #Value of first semi minor to semi major axis ratio
ca = 1   #Value of second semi minor to semi major axis ratio
N = 15   #Number of Particles in the simulation
dt = 10  #Time Step

#Initialize the particles and the velocities
particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)
velocities = np.zeros_like(particles) 

#Calculate the first velocity half step
densityField, _,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size=32)
phi = PoissonEquation.IsolatedMass(densityField)
velocities = dt/2 * IntegrationMethods.F(particles, IntegrationMethods.gradient(phi))



""" Update the next animation frame
Inputs:
 ------
 frame: value
     the number of times the update function needs to run

 Returns:
 -------
 scat,: Matplotlib.pyplot.plot.figure
     the plot showing particle position
"""
def update(frame):
    global particles, velocities
    particles, velocities = IntegrationMethods.xnext(particles, velocities, dt)
    scat._offsets3d = (particles[:, 0], particles[:, 1], particles[:, 2])
    return scat,

#Initialize the new plot
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
