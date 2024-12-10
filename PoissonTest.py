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

#central point Poisson test 
testpoint = np.array([0,0,0])


arr = particles #
#arr=np.array([testpoint])

densityField, x,y,z = InitializePoints.CreateDensityField(center,particles, grid_size)
InitializePoints.plotInitialPoints(arr)

#densityField, x,y,z = InitializePoints.CreateDensityField(center, arr, grid_size)
#testpotential=PoissonEquation.PoissonsEq(densityField)   

 
testpotential=PoissonEquation.IsolatedMass(densityField)   

#this creates the new axes for the extended mesh
x2 = np.linspace(-.5, 1.5, 64, endpoint=True) + center[0]
y2 = np.linspace(-.5, 1.5, 64, endpoint=True) + center[1]
z2 = np.linspace(-.5, 1.5, 64, endpoint=True) + center[2]

#This plots the potential using the density plotter function, so the axis label is wrong
InitializePoints.PlotTestFields(testpotential, x2, y2, z2, 32)






''' 
particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)
densityField2, x,y,z = InitializePoints.CreateDensityField(center, a, ba, ca, particles, grid_size)
InitializePoints.plotInitialPoints(arr)
testpotentialgauss=PoissonEquation.PoissonsEq(densityField)   
PoissonEquation.PlotPotential2D(testpotentialgauss, x,y,z,'x',int(grid_size/2))
PoissonEquation.PlotPotential2D(testpotentialgauss, x,y,z,'y', int(grid_size/2))
PoissonEquation.PlotPotential2D(testpotentialgauss, x,y,z,'z', int(grid_size/2))
'''
