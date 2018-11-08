# Title: Imitation Exponential Decay Gradient
# Author: Austin Foster
# Date: November 7, 2018
# Description:  This code produces a thermal gradient data set that
#               matches the exponential decay that Dr. Tolley's code
#               uses.

# Import modules
import numpy as np
import matplotlib.pyplot as plt

# Geometry Parameters
plateLen = 10.0
colLength = 1000.0 # measured in cm
nCorners = 200 # number of switchbacks
TempGranularity = 1.0e2 # data points/cm

# Create full precision x array
x_Full = np.linspace(0,colLength,colLength*TempGranularity)

# Define the exponenetial decay
decayConst = 0.5
TRamp = 50.0
expoDecay = np.exp(-np.arange(1, (plateLen*TempGranularity)+1)/\
                  TempGranularity*decayConst)*TRamp

# Compute column corner geometry
lenCross = (colLength-plateLen)/nCorners  #This gives the length of each
                                            # isotherm.
seqStart = 1
cornerSpacing = 1.0*(len(expoDecay)-seqStart)/(nCorners)    # This gives the 
corner = np.round(np.arange(seqStart,len(expoDecay)+1, cornerSpacing)) #indexes the position of each of the corners
            # on the exponential decay curve.
print(corner)
TempVec_Full = np.array([])
for jj in range(0,nCorners):
    #This gives the temperature at each of the corners
    TempVec_Full = np.append(TempVec_Full,expoDecay[int(corner[jj]-1):int(corner[jj+1])])
    #  This determines the temperature at each point along the isotherm.
    TempVec_Full = np.append(TempVec_Full, np.repeat(expoDecay[int(corner[jj+1]-1)],lenCross*TempGranularity))
TempVec_Full = TempVec_Full[0:int(colLength*TempGranularity)]

# ---- Low Resolution Gradient ---- #
TempVec_Lo = np.array([])
for jj in range(0,nCorners):
    #This gives the temperature at each of the corners
    TempVec_Lo = np.append(TempVec_Lo,[expoDecay[int(corner[jj]-1)],expoDecay[int(corner[jj])]])
    #  This determines the temperature at each point along the isotherm.
##    TempVec_Lo = np.append(TempVec_Lo, np.repeat(expoDecay[int(corner[jj+1]-1)],2))
plt.plot(TempVec_Lo,'r.-')
plt.show()
##TempVec_Lo = TempVec_Lo[0:int(colLength*TempGranularity)]

# Compare gradients
plt.plot(x_Full,TempVec_Full,'k.-')
plt.show()
