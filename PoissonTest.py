import InitializePoints
import PoissonEquation
import numpy as np




center = [0,0,0]
a = 4
ba = 1
ca = 1
N = 32**3
grid_size = 32



a = .2
ba = 1
ca = 1
N = 32**3
    
particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)

#central point discrete Poisson eq test (no Green's function) 
def central_potential_test():
 center = [0,0,0]
 grid_size = 32
    
 x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
 y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
 z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]
    
#Place one particle at exactly hte middle point of the grid
 particles = np.array([[x[16],y[16],z[16]]])
 densityField,_,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)
 potential=PoissonEquation.discretepoisson(densityField)
 green=PoissonEquation.IsolatedMass(densityField)
#analytic expression for Poisson's eq:
 analytic= np.zeros((32,32,32))
 for i in range(0,32):
  for j in range(0,32):
   for k in range(0,32):
    r2=((16-i)**2 +(16-j)**2 +(16-k)**2)
    if r2 ==0:
     analytic[i][j][k]=-32
    else:
     analytic[i][j][k]=-32/np.sqrt(r2)
 print(potential)
 print(analytic)    
 print(np.sum(potential-analytic) /32**3) 
 InitializePoints.PlotTestFields(densityField, x,y,z, grid_size)
# InitializePoints.PlotTestFields(potential, x,y,z, grid_size)
# InitializePoints.PlotTestFields(green, x,y,z, grid_size)
 PoissonEquation.PlotTestPotential(analytic, x,y,z, grid_size)
 PoissonEquation.PlotTestPotential(potential, x,y,z, grid_size)
 PoissonEquation.PlotTestPotential(green, x,y,z, grid_size)
 
# assert (not all.(analytic==potential) ) 
  
central_potential_test() 

     
#central point discrete Poisson eq test using Green's function convolution
def centralpoint_potential_green_function_test():
 center = [0,0,0]
 grid_size = 32
    
 x = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[0]
 y = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[1]
 z = np.linspace(-0.5, 0.5, grid_size, endpoint=True) + center[2]

#Place one particle at exactly hte middle point of the grid
 particles = np.array([[x[16],y[16],z[16]]])
 densityField,_,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)
 potential=PoissonEquation.discretepoisson(densityField)
    
    
   
    
#Place one particle at exactly the middle point of the grid
 particles = np.array([[x[16],y[16],z[16]]])
 densityField,_,_,_ = InitializePoints.CreateDensityField(center, particles, grid_size)
 potential=PoissonEquation.IsolatedMass(densityField) 

def spherical_potential_green_function_test():
 hi=1
#Potential from two oint particles:
