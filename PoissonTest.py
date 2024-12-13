import InitializePoints
import PoissonEquation
import numpy as np


    
particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)

#central point discrete Poisson eq test 
def central_potential_test():
 center = [0,0,0]
 grid_size = 32

#define grid axis values    
 x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
 y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
 z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]
    
#Place one particle at exactly the middle point of the grid, obtain
 particles = np.array([[x[16],y[16],z[16]]])
 densityField,_,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)
 potential=PoissonEquation.discretepoisson(densityField)
 green=PoissonEquation.IsolatedMass(densityField)
 
#define analytic expression for Poisson's eq:
 analytic= np.zeros((32,32,32))
 for i in range(0,32):
  for j in range(0,32):
   for k in range(0,32):
    r2=((16-i)**2 +(16-j)**2 +(16-k)**2)
    if r2 ==0:
     analytic[i][j][k]=-32
    else:
     analytic[i][j][k]=-32/np.sqrt(r2)

 #InitializePoints.PlotTestFields(densityField, x,y,z, grid_size)
# InitializePoints.PlotTestFields(potential, x,y,z, grid_size)
# InitializePoints.PlotTestFields(green, x,y,z, grid_size)
 #PoissonEquation.PlotTestPotential(analytic, x,y,z, grid_size)
 #PoissonEquation.PlotTestPotential(potential, x,y,z, grid_size)
 #PoissonEquation.PlotTestPotential(green, x,y,z, grid_size)
 
 #Test that analytic solution and green's function convolution are nearly identical
 assert not (np.abs(analytic-potential)<1e-7).all()
 
 assert (np.abs(analytic- green)<1e-7).all()


#Spherical potential test
def spherical_potential_green_function_test():
 center = [0,0,0]
 grid_size = 32
   
#define grid axis values     
 x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
 y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
 z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]
    
#create density field and potential from Green's function     
 densityField=InitializePoints.CreateSphericalShellDensityField(grid_size = 32, radius_sq = 0.3**2, tol = 0.02)
 green=PoissonEquation.IsolatedMass(densityField)
 
#define analytic expression for potential from density field:
 analytic= np.zeros((32,32,32))
 #obtain particle indices
 densidx = np.nonzero(densityField)
 densidx2 = np.zeros(densidx[0].shape)
 particles2=[]
 for num in range(len(densidx2)):
  part =[]
  for num2 in range(len(densidx)):
   part.append(densidx[num2][num])
  particles2.append(part) 
 partarr=np.array(particles2)
# print(partarr)
 #set potential values by finding the potential from each point source and superposing them
 for particle in partarr:
  particle_pot= np.zeros((32,32,32))
  for i in range(0,32):
   for j in range(0,32):
    for k in range(0,32):
     r2=((i-particle[0])**2 +(j-particle[1])**2 +(k-particle[2])**2)
     #   r2=((dens[0]-x[i])**2 +(dens[1]-y[j])**2 +(dens[2]-z[k])**2)
     if r2 ==0:
      particle_pot[i][j][k]=-32
     else:
      particle_pot[i][j][k]=-32/np.sqrt(r2)
  analytic += particle_pot     


 #InitializePoints.PlotTestFields(densityField, x,y,z, grid_size)
# InitializePoints.PlotTestFields(potential, x,y,z, grid_size)
# InitializePoints.PlotTestFields(green, x,y,z, grid_size)
 #PoissonEquation.PlotTestPotential(analytic, x,y,z, grid_size)
 ##PoissonEquation.PlotTestPotential(potential, x,y,z, grid_size)
 PoissonEquation.PlotTestPotential(green, x,y,z, grid_size)
 
# print(green)
# print(analytic) #
 print(np.abs(analytic- green)) 
 
 #Test that analytic solution and green's function convolution are nearly identical
 assert (np.abs(analytic - green)<1e-7).all()

 
 
#Potential from two point particles
def two_part_potential_green_function_test():
 center = [0,0,0]
 grid_size = 32

#define grid axis values    
 x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
 y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
 z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]
    
#Place one particle at exactly the middle point of the grid
 particles = np.array([[x[16],y[2],z[16]],[x[16],y[29],z[16]]])
 densityField,_,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)
 potential=PoissonEquation.discretepoisson(densityField)
 green=PoissonEquation.IsolatedMass(densityField)
#defne analytic expression for Poisson's eq:
 analytic= np.zeros((32,32,32))
 #obtain particle indices 
 densidx = np.nonzero(densityField)
 densidx2 = np.zeros(densidx[0].shape)
 particles2=[]
 
 for num in range(len(densidx2)):
  part =[]
  for num2 in range(len(densidx)):
   part.append(densidx[num2][num])
  particles2.append(part) 
 partarr=np.array(particles2)
 #print(partarr)
 #set potential values by finding the potential from each point source and superposing them 
 for particle in partarr:
  particle_pot= np.zeros((32,32,32))
  for i in range(0,32):
   for j in range(0,32):
    for k in range(0,32):
     r2=((i-particle[0])**2 +(j-particle[1])**2 +(k-particle[2])**2)
     #   r2=((dens[0]-x[i])**2 +(dens[1]-y[j])**2 +(dens[2]-z[k])**2)
     if r2 ==0:
      particle_pot[i][j][k]=-32
     else:
      particle_pot[i][j][k]=-32/np.sqrt(r2)
  analytic += particle_pot     


 #InitializePoints.PlotTestFields(densityField, x,y,z, grid_size)
# InitializePoints.PlotTestFields(potential, x,y,z, grid_size)
# InitializePoints.PlotTestFields(green, x,y,z, grid_size)
 PoissonEquation.PlotTestPotential(analytic, x,y,z, grid_size)
 #PoissonEquation.PlotTestPotential(potential, x,y,z, grid_size)
 PoissonEquation.PlotTestPotential(green, x,y,z, grid_size)
 
 #print(green)
 #print(analytic) #
 #print(np.abs(analytic- green)) 
 
 #Test that analytic solution and green's function convolution are nearly identical
 assert (np.abs(analytic - green)<1e-7).all()

# assert (not all.(analytic==potential) ) 


central_potential_test() 
spherical_potential_green_function_test()
two_part_potential_green_function_test()




