import InitializePoints
import PoissonEquation
import numpy as np


center = [0,0,0]
a = 4
ba = 1
ca = 1
N = 32**3
grid_size = 32


#central point Poisson test 
testpoint = np.array( [0, 0, 0])
arr = np.array([testpoint])

densityField, x,y,z = InitializePoints.CreateDensityField(center, a, ba, ca, arr, grid_size)
InitializePoints.plotInitialPoints(arr)
densityField, x,y,z = InitializePoints.CreateDensityField(center, a, ba, ca, arr, grid_size)
testpotential=PoissonEquation.PoissonsEq(densityField)   
PoissonEquation.PlotPotential2D(testpotential, x,y,z,'x',int(grid_size/2))
PoissonEquation.PlotPotential2D(testpotential, x,y,z,'y', int(grid_size/2))
PoissonEquation.PlotPotential2D(testpotential, x,y,z,'z', int(grid_size/2))


''' 
particles = InitializePoints.initializeGaussianPoints(center, a, ba, ca, N)
densityField2, x,y,z = InitializePoints.CreateDensityField(center, a, ba, ca, particles, grid_size)
InitializePoints.plotInitialPoints(arr)
testpotentialgauss=PoissonEquation.PoissonsEq(densityField)   
PoissonEquation.PlotPotential2D(testpotentialgauss, x,y,z,'x',int(grid_size/2))
PoissonEquation.PlotPotential2D(testpotentialgauss, x,y,z,'y', int(grid_size/2))
PoissonEquation.PlotPotential2D(testpotentialgauss, x,y,z,'z', int(grid_size/2))
'''
