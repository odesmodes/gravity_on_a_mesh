#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 20:22:55 2024

@author: vedhasyamuvva
"""

import numpy as np
import InitializePoints
import PoissonEquation

def gradient(phi, spacing = 1/32):
    grad = np.gradient(phi, spacing)
    return np.array([grad[0],grad[1],grad[2]])

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


def xnext(x, v, dt): 
    #RECALCULATE PHI
    densityField, _,_,_ = InitializePoints.CreateDensityField([0,0,0], x, grid_size = 32)
    phi = PoissonEquation.IsolatedMass(densityField)
    gradArr = gradient(phi)
    
    # Calculate gradArr with new distribution

    # VERLET METHOD FROM WIKIPEDIA
    
    F_old = F(x, gradArr)
    x_new = x + v * dt + 0.5 * F_old * dt**2
    gradArr_new = gradient(phi)
    F_new = F(x_new, gradArr_new)
    v_new = v + 0.5 * (F_old + F_new) * dt
    
    
    #Verlet METHOD from CLASS NOTES
    """
    v_new = v + F(x, gradArr) * dt/2
    #print("x: ", np.shape(x), np.shape(v_new))
    x_new = x + v_new * dt
    """
    return x_new, v_new

def NetAngularMomentum(particles):
    R_z = np.array([[0, -1, 0],
                    [1,  0, 0],
                    [0,  0, 1]])
    velocities = particles @ R_z.T
    
    magnitudes = np.linalg.norm(velocities, axis=1, keepdims=True)
    magnitudes[magnitudes == 0] = 1
    
    return velocities / magnitudes

def MassWithinRadius(densityField,r, tol = 1e-6):
    grid_size = 32
    center = [0,0,0]

    x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
    y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
    z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]
    
    # Generate grid coordinates
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    # Calculate squared distance from the center
    distance_sq = X**2 + Y**2 + Z**2

    # Identify cells within the specified radius
    within_radius = (distance_sq <= r**2 + tol)

    # Calculate total mass (density * volume) for cells within the radius
    total_mass = np.sum(densityField[within_radius])

    return total_mass

def NoEvolution(particles):
    densityField, _,_,_ = InitializePoints.CreateDensityField([0,0,0], particles , grid_size = 32)
    
    distance = np.linalg.norm(particles, axis=1, keepdims=True)
    refVec = np.array([0,0,1])
    refVec2 = np.array([1,0,0])
    velocities = np.empty_like(particles)
    
    for i in range (np.shape(particles)[0]):
        ptcl = particles[i]
        vhat = np.cross(ptcl, refVec)/distance[0]
        if np.allclose(vhat, np.array([0, 0, 0])): vhat = np.cross(ptcl, refVec2)
        velocities[i] = MassWithinRadius(densityField,distance[i]) / distance[i] * vhat
    
    return velocities
    
    
    
    
    