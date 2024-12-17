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
    """
    Calculates the gradient of a scalar field (phi) in 3D using finite differences.

    Inputs:
    ------
    phi : NumPy array
        The scalar field for which the gradient is calculated.

    spacing : float, optional
        The spacing between grid points (default is 1/32).

    Returns:
    -------
    grad : NumPy array
        The gradient of the scalar field phi, represented as a 3D vector field.
    """
    grad = np.gradient(phi, spacing)
    return np.array([grad[0], grad[1], grad[2]])


def F(particles, gradArr, spacing=1/32):
    """
    Calculates the forces on particles based on the gradient of a potential field.

    Inputs:
    ------
    particles : NumPy array
        Array of particle positions.

    gradArr : NumPy array
        Array representing the gradient field at the particle positions.

    spacing : float, optional
        The spacing between grid points (default is 1/32).

    Returns:
    -------
    forces : NumPy array
        The forces acting on each particle.
    """
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
    """
    Calculates the next position and velocity of particles using the Verlet integration method.

    Inputs:
    ------
    x : NumPy array
        Current positions of the particles.

    v : NumPy array
        Current velocities of the particles.

    dt : float
        Time step for the integration.

    Returns:
    -------
    x_new : NumPy array
        Updated positions of the particles.

    v_new : NumPy array
        Updated velocities of the particles.
    """
    # RE-CALCULATE PHI
    densityField, _, _, _ = InitializePoints.CreateDensityField([0, 0, 0], x, grid_size=32)
    phi = PoissonEquation.IsolatedMass(densityField)
    gradArr = gradient(phi)
    
    # Calculate gradArr with new distribution

    # Verlet method from Wikipedia
    F_old = F(x, gradArr)
    x_new = x + v * dt + 0.5 * F_old * dt**2
    gradArr_new = gradient(phi)
    F_new = F(x_new, gradArr_new)
    v_new = v + 0.5 * (F_old + F_new) * dt
    
    return x_new, v_new


def NetAngularMomentum(particles):
    """
    Calculates the net angular momentum of particles in the system.

    Inputs:
    ------
    particles : NumPy array
        Array of particle positions and velocities.

    Returns:
    -------
    velocities : NumPy array
        The angular velocities of the particles.
    """
    R_z = np.array([[0, -1, 0],
                    [1,  0, 0],
                    [0,  0, 1]])
    velocities = particles @ R_z.T
    
    magnitudes = np.linalg.norm(velocities, axis=1, keepdims=True)
    magnitudes[magnitudes == 0] = 1
    
    return velocities / magnitudes


def MassWithinRadius(densityField, r, tol=1e-6):
    """
    Calculates the total mass within a given radius from the center based on the density field.

    Inputs:
    ------
    densityField : NumPy array
        The density field of the system.

    r : float
        The radius within which to calculate the mass.

    tol : float, optional
        A small tolerance value to account for the thickness of the shell (default is 1e-6).

    Returns:
    -------
    total_mass : float
        The total mass enclosed within the radius r.
    """
    grid_size = 32
    center = [0, 0, 0]

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
    h = 1/grid_size
    total_mass = np.sum(densityField[within_radius]) * h**3

    return total_mass


def NoEvolution(particles):
    """
    Calculates the velocities of particles assuming no evolution in the system, based on their positions.

    Inputs:
    ------
    particles : NumPy array
        Array of particle positions.

    Returns:
    -------
    velocities : NumPy array
        The calculated velocities of the particles assuming no evolution.
    """
    densityField, _, _, _ = InitializePoints.CreateDensityField([0, 0, 0], particles, grid_size=32)
    
    distance = np.linalg.norm(particles, axis=1, keepdims=True)
    refVec = np.array([0, 0, 1])
    refVec2 = np.array([1, 0, 0])
    velocities = np.empty_like(particles)
    
    for i in range(np.shape(particles)[0]):
        if (distance[i] == 0):
            velocities[i] = np.array([0, 0, 0])
        else:
            ptcl = particles[i]
            vhat = np.cross(ptcl, refVec) / distance[i]
            if np.allclose(vhat, np.array([0, 0, 0])): 
                vhat = np.cross(ptcl, refVec2)
            velocities[i] = MassWithinRadius(densityField, distance[i]) / distance[i] * vhat
    
    return velocities
