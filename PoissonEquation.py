import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#Section 3

def PoissonsEq(density):
 '''Takes density as an input, solves Poissons eq using Fourier Transformed density, inverse Fourier Transforms result to return solution for potential'''
 density_dft = sp.fft.fftn(density)
 potential_dft = np.empty((33,33,33),dtype = 'complex_')
 for x in range(1,len(density_dft)):
  for y in range(1,len(density_dft[x])):
   for z in range(1,len(density_dft[x][y])):
    potential_dft[x][y][z]=1/(2*(np.cos(2*x*np.pi/len(density_dft))+np.cos(2*y*np.pi/len(density_dft[x]))+np.cos(2*z*np.pi/len(density_dft[x][y])) -3))*4*np.pi*density_dft[x][y][z]
#  print('potential dft: '+str(potential_dft))
 potential_ift = sp.fft.ifftn(potential_dft)
 potential=potential_ift.real
 print('potential: '+str(potential))

 return potential
 
def PlotPotential2D(potential, x,y,z, axis, value):
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


'''

def IsolatedMass(density):
"""Solves Poisson's equation for an isolated mass distribution using a convolution of the Fourier transforms of the mass m"""
#create extended mesh to isolate mass 
 iso_den= density.resize((66,66,66)
#create Green's function with restricted range
 for x in range(1,33):
  for y in range(1,33):
   for z in range(1,33):
    G[x][y][z]=1/np.sqrt(x**2 + y**2+ z**2)
 
    
  G[0][0][0]=1
 
#FT Green's function
 G_dft= sp.fft.fftn
 
#
#FT density
 density_ft=sp.fft.fftn(density) 
  
 #
 return potential
 
 

 '''
 
