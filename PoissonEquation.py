import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#Section 3

def PoissonsEq(density):
 '''Takes density as an input, solves Poissons eq using Fourier Transformed density, inverse Fourier Transforms result to return solution for potential'''
#Take fft of density 
 density_dft = sp.fft.fftn(density)

#Define array to hold F-transformed potential  
 potential_dft = np.empty((32,32,32),dtype = 'complex_')
 
#loop over density FFT and multiply by expression discrete transform  
 for x in range(0,len(density_dft)):
  for y in range(0,len(density_dft[x])):
   for z in range(0,len(density_dft[x][y])):
    if 2*(np.cos(2*x*np.pi/len(density_dft))+np.cos(2*y*np.pi/len(density_dft[x]))+np.cos(2*z*np.pi/len(density_dft[x][y])) -3)==0:
    #if expression is 0, set potential dft to 0
     potential_dft[x][y][z]=0
    else: 
     potential_dft[x][y][z]=1/(2*(np.cos(2*x*np.pi/len(density_dft))+np.cos(2*y*np.pi/len(density_dft[x]))+np.cos(2*z*np.pi/len(density_dft[x][y])) -3))*4*np.pi*density_dft[x][y][z]
    
#take inverse FT of potental FT
 potential_ift = sp.fft.ifftn(potential_dft)
 
#take real part of potential and return 
 potential=potential_ift.real
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
 iso_den=np.zeros((64,64,64))

#create array for Green's function with restricted range
 G=np.zeros((64,64,64))

#loop over points and set to 1/|r| where r is the distance from each outermost vertex of the extended cubic mesh, within each of the 8 smaller cubic regions containing these vertices. The active region in which the mass distribution is contained is x,y,z <32
 for x in range(0,64):
  for y in range(0,64):
   for z in range(0,64):
    if x<32:
     xval=x
    else:
     xval=(64-x) 
    if y<32:
     yval=y
    else:
     yval=(64-y)   
    if z<32:
     zval =z
    else:
     zval=(64-z)
    r2=(xval**2 +yval**2+zval**2)
    if r2==0:
     G[x][y][z]=1
    else: 
     G[x][y][z]=1/np.sqrt(r2)
#set the iso_den to density value with in active region     
    if x<32 and y < 32 and z < 32: 
     iso_den[x][y][z]=density[x][y][z]   


#FT Green's function
 G_ft= sp.fft.fftn(G)
#FT density
 density_ft=sp.fft.fftn(iso_den) 

#Define Potential FT array 
 potential_ft=np.zeros((64,64,64),dtype = 'complex_')
 
#set potential ft to G_ft *density_ft 
 for x in range(0,64):
  for y in range(0,64):
   for z in range(0,64):
    potential_ft[x][y][z] = G_ft[x][y][z]*density_ft[x][y][z]

#invert potential FT and take real part  
 potential_ift = sp.fft.ifftn(potential_ft)
 potential=potential_ift.real

 return potential



