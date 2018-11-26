import numpy as np
import math
colors  = ['red','green','blue','grey30']
markers = np.array(["o","^","P","x","D","v","s"])
j = 2                       # Number of compounds
nMol = 100                  # Number of molecules per compound used in simulation
sigmaNot = 10000*math.sqrt(1.0/12)
gamma1 = 5.113e-3           # Molecular diffusion coefficient when using both temp and pressure.
taylor = 1.647e-1           # Dispersion coefficient due to resistence to transport and taylor dispersion.
gamma3 = 7.676e-8           # Resistence to flow for adsorption/deporption (no pressure variable)
gamma2 = 1.217e-9           # Golay's formula for resistence to transport when using both temp and pressure in diffusion expression.
chkSm = 0
X = np.array([9.0,10.0])    # Alkane length
delta_t = 0.0001            # Time step
R = 8.3144621               # Boltzman's constant [Joules/ Kelvin mole]
colLength = 1000.0          # Column length [cm]
T0 = 290.0                  # Initial temperature [Kelvin, K]
diameter = 0.0100           # Inner diameter of the column [cm]
p_i = 57.16                 # Column inlet pressure by gage measured in psi.  For true inlet must add atmosphere pressure.
p_o = 13.8                  # This is the outlet pressure or Atmospheric pressure in psi.
w = 6.0                     # velocity of gradient [cm/sec]
tpRate = .25                # degrees C per second. Rate at which the column is heated after a gradient is established.
TempGranularity = 1.0e2     # granularity of temperatures across column. data points/cm
powerAce = 100.0