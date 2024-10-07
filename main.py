import InitializePoints

center = [0,0,0]
a = 1
ba = 2
ca = 3
N = 32**3


arr = InitializePoints.initializeGaussianPoints(center, a, ba, ca)

InitializePoints.plotInitialPoints(arr)
