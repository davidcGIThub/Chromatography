import numpy as np
import math

gradientFileName = "temperatureTest.csv"
positionFileName = "positionTest.csv"
timeFileName = "timeTest.csv"
posUnits = 1 # assume cm, (if meters then posUnits= 100, if mm then posUnits = .1 etc . .)

colors  = ['red','green','blue','grey30']
markers = np.array(["o","^","x","D","v","s"])
j = 2                     # Number of compounds
X = np.array([9.0,10.0])    # Alkane length
nMol = 100                  # Number of molecules per compound used in simulation
gamma1 = 5.113e-3           # Molecular diffusion coefficient when using both temp and pressure.
gamma3 = 7.676e-8           # Resistence to flow for adsorption/deporption (no pressure variable)
gamma2 = 1.217e-9           # Golay's formula for resistence to transport when using both temp and pressure in diffusion expression.
delta_t = 0.0001            # Time step
R = 8.3144621               # Boltzman's constant [Joules/ Kelvin mole]
diameter = 0.0100           # Inner diameter of the column [cm]
p_i = 57.16                 # Column inlet pressure by gage measured in psi.  For true inlet must add atmosphere pressure.
p_o = 13.8                  # This is the outlet pressure or Atmospheric pressure in psi.
tempGranularity = 1.0e2     # granularity of temperatures across column. data points/cm
powerAce = 100.0
injectionWidth = 2          #  Width of injection plug at time zero in cm.
