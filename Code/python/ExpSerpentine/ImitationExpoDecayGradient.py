# Title: Imitation Exponential Decay Gradient
# Author: Austin Foster
# Date: November 7, 2018
# Description:  This code produces a thermal gradient data set that
#               matches the exponential decay that Dr. Tolley's code
#               uses.

# Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

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

TempVec_Full = np.array([])
for jj in range(0,nCorners):
    #This gives the temperature at each of the corners
    TempVec_Full = np.append(TempVec_Full,expoDecay[int(corner[jj]-1):int(corner[jj+1])])
    #  This determines the temperature at each point along the isotherm.
    TempVec_Full = np.append(TempVec_Full, np.repeat(expoDecay[int(corner[jj+1]-1)],lenCross*TempGranularity))
TempVec_Full = TempVec_Full[0:int(colLength*TempGranularity)]

# ---- Low Resolution Gradient ---- #
TempVec_Lo = np.array([])
x_Lo = np.array([0.0])
flag = 1
for jj in range(0,nCorners-1):
    #This gives the temperature at each of the corners
    TempVec_Lo = np.append(TempVec_Lo,[expoDecay[int(corner[jj]-1)],expoDecay[int(corner[jj+1]-1)]])
    #  This determines the temperature at each point along the isotherm.
##    TempVec_Lo = np.append(TempVec_Lo, np.repeat(expoDecay[int(corner[jj+1]-1)],2))
    # Add x data
    x_Lo = np.append(x_Lo,x_Lo[-1]+(corner[jj+1]-corner[jj])/TempGranularity)
    if (corner[jj]%5==0 and flag):
        x_Lo = np.append(x_Lo,x_Lo[-1]+lenCross+2/TempGranularity)
        flag = 0
    else:
        x_Lo = np.append(x_Lo,x_Lo[-1]+lenCross+1/TempGranularity)
x_Lo = x_Lo[0:len(x_Lo)-1]
##plt.plot(x_Lo,TempVec_Lo,'r.-')
##TempVec_Lo = TempVec_Lo[0:int(colLength*TempGranularity)]

### Compare gradients
##plt.plot(x_Full,TempVec_Full,'k.')
##plt.show()

# Time increment parameters
runTime_n = 1500000
delta_t_Tolley = 0.0001
t_Tot = delta_t_Tolley*runTime_n
delta_t = 1

# Temperature program ramp rate parameters
T0 = 290.0
tpRate = 0.25    # Degrees C per second

# Create time array
time = np.arange(0,t_Tot+delta_t,delta_t)

# Open temperature matrix creation loop
TempMat = np.zeros((len(time),len(x_Lo)))
for i in range(len(time)):
    T_ramp = TempVec_Lo + i*delta_t*tpRate + T0
    TempMat[i,:] = T_ramp

# Create grid for plotting
xG,tG = np.meshgrid(x_Lo,time)

# Save data
save = 1
if save:
    np.savetxt('time.csv',time,delimiter=',')
    np.savetxt('colPosition.csv',x_Lo,delimiter=',')
    np.savetxt('temperature.csv',TempMat,delimiter=',')

# Display 3D plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_wireframe(xG,tG,TempMat)
ax.set_xlabel('Column Position (m)')
ax.set_ylabel('Time (s)')
ax.set_zlabel('Temperature (C)')
plt.show()
