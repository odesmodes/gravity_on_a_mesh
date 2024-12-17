import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#Section 3

def discretepoisson(density):
    ''' 
    Solves the Fourier-transformed version of Poisson's equation for a discrete mass distribution.

    Takes the density field as input, applies the Fourier transform, multiplies by an appropriate factor 
    (the discrete equivalent of k**-2), and then takes the inverse Fourier transform to compute the potential.

    Parameters:
    ----------
    density : numpy.ndarray
        The input 3D density field for which the potential is to be computed.

    Returns:
    -------
    potential : numpy.ndarray
        The 3D potential field corresponding to the input density.
    '''
    # Take fft of density 
    density_dft = sp.fft.fftn(density)

    # Define array to hold F-transformed potential  
    potential_dft = np.empty((32, 32, 32), dtype='complex_')

    # Loop over density FFT and multiply by expression for discrete transform  
    for x in range(0, len(density_dft)):
        for y in range(0, len(density_dft[x])):
            for z in range(0, len(density_dft[x][y])):
                # Multiply by factor (discrete equivalent of k**-2)
                if 2 * (np.cos(2 * x * np.pi / len(density_dft)) + np.cos(2 * y * np.pi / len(density_dft[x])) + np.cos(2 * z * np.pi / len(density_dft[x][y])) - 3) == 0:
                    # If expression is 0, set potential dft to 0
                    potential_dft[x][y][z] = 0
                else: 
                    potential_dft[x][y][z] = 1 / (32**2 * 2 * (np.cos(2 * x * np.pi / len(density_dft)) + np.cos(2 * y * np.pi / len(density_dft[x])) + np.cos(2 * z * np.pi / len(density_dft[x][y])) - 3)) * 4 * np.pi * density_dft[x][y][z]

    # Take inverse FT of potential FT
    potential_ift = sp.fft.ifftn(potential_dft)

    # Take real part of potential and return 
    potential = potential_ift.real
    return potential

 
def PlotPotential2D(potential, x, y, z, axis, value):
    ''' 
    Plots the 2D cross-sectional view of the potential at a specific slice along the specified axis.

    Parameters:
    ----------
    potential : numpy.ndarray
        The 3D potential field to be plotted.

    x, y, z : numpy.ndarray
        The 1D arrays representing the spatial coordinates along each axis.

    axis : str
        The axis along which to take the slice ('x', 'y', or 'z').

    value : int
        The index at which to slice the potential along the specified axis.
    '''
    if axis == 'y':
        extent = [x[0], x[-1], x[0], x[-1]]  # Define the physical coordinates for the plot
        plt.imshow(potential[:, value, :], origin='lower', cmap='inferno', extent=extent, aspect='auto')
        plt.title(f'Potential at Y={y[value]}')
        plt.ylabel('X-axis')
        plt.xlabel('Z-axis')
    elif axis == 'x':
        extent = [x[0], x[-1], x[0], x[-1]]
        plt.imshow(potential[value, :, :], origin='lower', cmap='inferno', extent=extent, aspect='auto')
        plt.title(f'Potential Slice at X={x[value]}')
        plt.ylabel('Y-axis')
        plt.xlabel('Z-axis')
    elif axis == 'z':
        extent = [x[0], x[-1], x[0], x[-1]]
        plt.imshow(potential[:, :, value], origin='lower', cmap='inferno', extent=extent, aspect='auto')
        plt.title(f'Potential Slice at Z={z[value]}')
        plt.ylabel('X-axis')
        plt.xlabel('Y-axis')

    plt.colorbar(label='Potential')
    plt.savefig(f"2DPotentialPlot{axis}_{value}.png")
    plt.show()

def PlotTestPotential(potential, x, y, z, grid_size):
    ''' 
    Plots the 2D cross-sectional views of the potential for all three axes (x, y, and z) at the middle slice.

    Parameters:
    ----------
    potential : numpy.ndarray
        The 3D potential field to be plotted.

    x, y, z : numpy.ndarray
        The 1D arrays representing the spatial coordinates along each axis.

    grid_size : int
        The size of the grid to calculate the center slice for plotting.
    '''
    PlotPotential2D(potential, x, y, z, "x", int(grid_size / 2))
    PlotPotential2D(potential, x, y, z, "y", int(grid_size / 2))
    PlotPotential2D(potential, x, y, z, "z", int(grid_size / 2))


def IsolatedMass(density):
    ''' 
    Solves Poisson's equation for an isolated mass distribution using the Fourier transform method with a Green's function.

    Takes the density field as input and computes the gravitational potential by performing Fourier transforms
    on both the density field and Green's function, then returning the resulting potential for the isolated mass.

    Parameters:
    ----------
    density : numpy.ndarray
        The 3D density field representing the mass distribution.

    Returns:
    -------
    potential : numpy.ndarray
        The 3D gravitational potential corresponding to the input density.
    '''
    # Create extended mesh to isolate mass 
    iso_den = np.zeros((64, 64, 64))

    # Create array for Green's function with restricted range
    G = np.zeros((64, 64, 64))

    # Loop over points and set to 1/|r| where r is the distance from each outermost vertex of the extended cubic mesh
    # within each of the 8 smaller cubic regions containing these vertices. The active region in which the mass distribution is contained is x,y,z < 32
    for x in range(0, 64):
        for y in range(0, 64):
            for z in range(0, 64):
                if x < 32:
                    xval = x / 32
                else:
                    xval = (64 - x) / 32 
                if y < 32:
                    yval = y / 32
                else:
                    yval = (64 - y) / 32  
                if z < 32:
                    zval = z / 32
                else:
                    zval = (64 - z) / 32
                r2 = (xval**2 + yval**2 + zval**2)
                if r2 == 0:
                    G[x][y][z] = 32
                else: 
                    G[x][y][z] = 1 / np.sqrt(r2)  
    # Set the iso_den to density value within active region     
    for x in range(0, 32):
        for y in range(0, 32):
            for z in range(0, 32): 
                iso_den[x][y][z] = density[x][y][z]   

    # FT Green's function
    G_ft = sp.fft.fftn(G)
    # FT density
    density_ft = sp.fft.fftn(iso_den) 

    # Define Potential FT array 
    potential_ft = np.zeros((64, 64, 64), dtype='complex')

    # Set potential ft to G_ft * density_ft 
    for x in range(0, 64):
        for y in range(0, 64):
            for z in range(0, 64):
                potential_ft[x][y][z] = G_ft[x][y][z] * density_ft[x][y][z]

    # Invert potential FT and take real part  
    potential_ift = sp.fft.ifftn(potential_ft)
    potential = potential_ift.real
    return -potential[0:32, 0:32, 0:32]
