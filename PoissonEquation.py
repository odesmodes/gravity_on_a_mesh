import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#Section 3

def PoissonsEq(density):
 '''Takes density as an input, solves Poissons eq using Fourier Transformed density, inverse Fourier Transforms result to return solution for potential'''
#Take fft of density 
 density_dft = sp.fft.fftn(density)
 potential_dft = np.empty((32,32,32),dtype = 'complex_')
 for x in range(1,len(density_dft)):
  for y in range(1,len(density_dft[x])):
   for z in range(1,len(density_dft[x][y])):
    potential_dft[x][y][z]=1/(2*(np.cos(2*x*np.pi/len(density_dft))+np.cos(2*y*np.pi/len(density_dft[x]))+np.cos(2*z*np.pi/len(density_dft[x][y])) -3))*4*np.pi*density_dft[x][y][z]
  print('potential dft: '+str(potential_dft))
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



def IsolatedMass(density):
 """Solves Poisson's equation for an isolated mass distribution using a convolution of the Fourier transforms of the mass m"""
#create extended mesh to isolate mass 
 iso_den= np.resize(density,(64,64,64))
#create Green's function with restricted range
 G=np.zeros((64,64,64))
#Greens' function on 32*3 grid 
 
 G32=np.zeros((32,32,32))
#loop over points except origin and set to 1/|r| 
 for x in list(range(1,16)) +list(range(17,32)):
  for y in list(range(1,16)) +list(range(17,32)):
   for z in list(range(1,16)) +list(range(17,32)):
    G[x][y][z]= 1/np.sqrt((16-x)**2 + (16-y)**2+(16-z)**2)   
#set origin to 1    
    G32[x][y][z]=G[x][y][z]
 G[16][16][16]=1
 G32[16][16][16]=1
 print(G32[16][16][16])
 print(G32[15][15][15])
 print(G32)
#FT Green's function
 G_dft= sp.fft.fftn(G)
 
#
#FT density
 density_ft=sp.fft.fftn(iso_den) 
 potential_dft = G_dft*density_ft
 potential_ift = sp.fft.ifftn(potential_dft)
 potential=potential_ift.real
 #
 return G32

 


