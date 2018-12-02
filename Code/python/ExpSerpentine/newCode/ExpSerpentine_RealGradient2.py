import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import Parameters_ExpSerpentine
import interpolateTemperature2 as interpT

########################
###   NEW EQUATION:  This is a step by step program to generate a random walk model
###     for a varietly of moving thermal gradient profiles.
###		This code was UPDATED in Riva with Samuel in May 2016
####################
TemperatureGradients = np.loadtxt('temperature.csv',delimiter = ',')
colPosition = np.loadtxt('colPosition.csv',delimiter = ',')
timeData = np.loadtxt('time.csv',delimiter = ',')
GData = interpT.GradientData(colPosition,timeData,TemperatureGradients)
colLength = colPosition[len(colPosition)-1] - colPosition[0]  # Column length [cm]
#####Parameters
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
diameter = 0.0100           # Inner diameter of the column [cm]
p_i = 57.16                 # Column inlet pressure by gage measured in psi.  For true inlet must add atmosphere pressure.
p_o = 13.8                  # This is the outlet pressure or Atmospheric pressure in psi.
w = 6.0                     # velocity of gradient [cm/sec]
tpRate = .25                # degrees C per second. Rate at which the column is heated after a gradient is established.
TempGranularity = 1.0e2     # granularity of temperatures across column. data points/cm
powerAce = 100.0
#####################

C = -7.22173 - 0.47406*X                # Formula for calculating C where X is the number of carbons.  Includes phase ratio
C = np.array([-11.1277,-11.2829])       # ????? Stick with this assignment ????
h = (2284.0 + 819*X) * 4.18442979       #  Calculation of enthalpy (Joules)
h = np.array([34926.27 , 35245.71])     # ????? stick with this assigment ?????

p_i = p_i*6894.757          # This converts pressure from psi units to Pascals
p_o = p_o*6894.757          # This converts pressure from psi units to Pascals
p_i = p_i + p_o             #This is the true inlet pressure.

Temp2PowVec = np.arange(0,1000+1/powerAce,1/powerAce)**1.646        #address the temperature /100, get the temp raised to the power
#################
#
#	The vector 'Temp2Pow' is sequence of temperatures raised to the 1.646 power
#		to adjust for the viscosity of helium.
#
#################
runTime = int( (timeData[len(timeData)-1] - timeData[0]) / delta_t)  #Number of potential steps to get compounds out of column.
pauseCount = int(1/delta_t) #number of frames to skip in plotting results.  The 'big.matrix'
                            # takes a snapshot at every 'pause.count' time for illustrating the
                            # progress of the molecules.
##### Start Standard Execution
x = np.zeros((j,nMol))
detector = np.zeros((j,nMol))         # This matrix keeps track of x-axis location of each
                                    # molecule of each analyte. Initially each is at point
                                    # '0' as indicated here.  The y-axis location is simply the
                                    # scaled value of the molecule's index number in the matrix x.
injectionWidth = 2  #  Width of injection plug at time zero in cm.
for i in range(0,j):
  x[i,] = np.random.choice(nMol, nMol, replace=False)*1.0*injectionWidth/nMol
# The above loop returns an order from 1 to n.mol for each molecule.  Since it is divided
# by n.mol, in essence we get each molecule uniquely with a value between 0 and 1.
# Multiplying by injection.width gives a plug that has the n.mol molecules uniformly distributed
# from 0 cm to injection.width cm at the beginning of the column.
#
#   Note that the standard deviation below is calculated from s^2 and not from the var
#   of the uniform. HDT
#
#   Remember, x[i,] is the vector of locations in x-axis of i-th analyte.  This location is used
#	at each step for each molecule to determine temperature and velocity.
#	This is also used to plot location of the molecules.  The y-axis values for plot are determined
#	from the matrix column index.  (Rows indicate analyte and column indicate individual molecules.)
currentTime = timeData[0]
resEnd = np.zeros((j,2))
eluted1 = 1
eluted2 = 1
detected1 = np.array([0])
detected2 = np.array([0])
##########  Temporary pseudo-do-loop for debugging.
# n <- 5
##########  This is where the do-loop for n starts.
totalTime1 = 0
totalTime2 = 0
spread1 = 0
spread2 = 0

peak_width = np.array([])
plot_res = np.array([])
#added_dc#############################################################
vel1_x = np.array([])
vel2_x = np.array([])
pos1_x = np.array([])
pos2_x = np.array([])
std1 = np.array([])
std2 = np.array([])
time = np.array([])

for n in range(0,runTime+1):
#for n in range(0,2):
    TempPlus = GData.getTemperatures(n*delta_t,TempGranularity)
    Tmax = np.amax(TempPlus)
    Tmin = np.amin(TempPlus)
    Tdelta = Tmax-Tmin
    temp = abs(x*TempGranularity+1)
    temp[temp >= len(TempPlus)] = len(TempPlus) - 1
    T_all = TempPlus[temp.astype(int)].flatten()
    C_all = C
    tmp_C_all = np.repeat([C_all],len(T_all)/len(C_all),axis=0).flatten('F')
    tmp_h = np.repeat([h],len(T_all)/len(h),axis = 0).flatten('F')
    k_all = np.exp(tmp_h/(R*T_all)+tmp_C_all)

       # Here I adjust the diffusion calculation according to the position of each molecule and the temperature
       #   at that location, taking adjusting the velocity resulting from the different pressure and viscosity values.  (HDT)
       # 	To do this there are several steps.
       #	First step is to determine a numerical value for the integral of temperature times viscosity.
    TSumPow = np.cumsum(Temp2PowVec[(TempPlus*powerAce).astype(int)-1]) / TempGranularity
       # integral of Temp^1.646
         # Note baseline viscosity cancels out.
         # T.sum.pow[length(T.sum.pow)]
         # is sum or integral from 0 to end.
    totInt = TSumPow[int(colLength*TempGranularity)-1]  #This is the integral from 0 to L of T^1.646.
    HeFactor = 0.000018662/(273.15**.646)  #Mult factor viscosity of He solvent.  This is delta.not in notes.
    CStar = (p_i**2-p_o**2)/(2*HeFactor*TSumPow[int(colLength*TempGranularity)-1])
         ## C.star is =to (pi^2-po^2)/(2*integral over L of T^1.646)
         #	I am doing the integral as simply the sum of the histogram pieces.  I could
         #	improve this step by using Simpson's rule or the Newton binomial trick if the
         #	Temperature can be modeled as a simple function of two variates. (HDT)
    velocity_x = ( (p_i**2-p_o**2)*TempPlus[temp.astype(int)] * (diameter/2)**2 / ((16*HeFactor*totInt )*np.sqrt(abs(p_i**2-(p_i**2-p_o**2)*TSumPow[temp.astype(int)]/totInt ))) )
       # Velocity updated by HDT on 31 March 2015
    p_x = np.sqrt(abs((p_i**2-(p_i**2-p_o**2)*TSumPow[temp.astype(int)]/totInt)) ).flatten()
       #This is the pressure at x
       #
       # Pressure updated by HDT on 31 March 2015
    TempPrssRat = (T_all**1.75)/p_x
    Dt = 2*(( TempPrssRat*gamma1 + ((1+6*k_all+11*k_all**2)/((1+k_all)**2))*velocity_x.flatten()**2*gamma2/(TempPrssRat) + gamma3*velocity_x.flatten()**2*k_all/(1+k_all) )/(1+k_all))
    sigmaDelta = np.sqrt(abs(Dt*delta_t))
    x = np.reshape((x.flatten() + ( velocity_x.flatten()*delta_t/(1 + k_all) + np.random.normal(0,sigmaDelta,j*nMol))),np.shape(x))
       # For the Milshtein correction term we have
       #  W.Lang <- rnorm(j*n.mol, mean = 0, sd = sigma.delta )
       # x <- x + ( velocity.x*delta.t +W.Lang +sigma.delta^2*((W.lang-delta.t)^2)/2 )/(1 + k.all)
    # if n%500==0:
    #     if np.max(x[0,]) > colLength:
    #         allDetected1 = x[0,][x[0,] > colLength]
    #         newDetected1 = allDetected1[allDetected1 != detected1]
    #         detected1 = np.concatenate(detected1, newDetected1)
    #         detector[0, newDetected1] = n*delta_t
    #     if np.max(x[1,]) > colLength:
    #         allDetected2 = x[1,][x[1,] > colLength]
    #         newDetected2 = allDetected2[allDetected2 != detected2]
    #         detected2 = np.concatenate(detected2, newDetected2)
    #         detector[1, newDetected2] = n*delta_t
    #         if len(detected_2)+len(detected_1)-2 >= 2*nMol:
    #             #break
    #             break
    if n % pauseCount == 0:
        #added DC ################################################
        vel1_x = np.append(vel1_x , np.mean(velocity_x[0]))
        vel2_x = np.append(vel2_x , np.mean(velocity_x[1]))
        pos1_x = np.append(pos1_x , np.mean(x[0]))
        pos2_x = np.append(pos2_x , np.mean(x[1]))
        standev1 = np.std(x[0,],axis=0,ddof=1)
        standev2 = np.std(x[1,],axis=0,ddof=1)
        std1 = np.append(std1 , standev1)
        std2 = np.append(std2 , standev2)

        if np.min(x) > colLength:
           #break
           break
        top_pks = np.zeros(j)
        pk_head = np.zeros(j)
        pk_tail = np.zeros(j)
        peakx = np.zeros(1)

        f,(ax1,ax2) = plt.subplots(2, sharex=True)
        plt.xlim(0,colLength)
        plt.ylim(Tmin-Tmax)
        ax1.yaxis.set_label_position("left")
        ax2.yaxis.set_label_position("left")
        ax1.set_ylabel("Temperature (K)")
        for i in range(0,j):
           lab = round(4*np.std(x[i,],axis=0,ddof=1) , 2)
           ax1.scatter(x[i,], 1.0*np.arange(1,nMol+1)/nMol*Tdelta+Tmin,color = colors[i],s = 0.5,  marker = markers[i],label =lab)
        ax1.legend()
        tmp = (np.std(x[0,],axis=0,ddof=1)+np.std(x[1,],axis=0,ddof=1))/2.0
        peak_width = np.append(peak_width, 4*tmp)
        res = (np.mean(x[0,])-np.mean(x[1,]))/(4.0*tmp)
        plot_res = np.append(plot_res, res)
        xx = np.linspace(0,colLength-1,(colLength-1)/.1 + 1)
        yy = (TempPlus[(xx*TempGranularity).astype(int)])
        ax1.plot(xx,yy,color = 'red')
        #abline(h=( (delta.t*mm*pause.count*b) )/200,col='red') #programmed temperature
        ax2.set_ylabel('Density')
        for i in range(0,j):
           sns.kdeplot(x[i,],color=colors[i],ax=ax2) #density plot
        currentTime = round(n*delta_t,3)
        intSec = int(currentTime)
        intDecimal = int((currentTime - intSec)*100)
        print(currentTime)
        time = np.append(time,currentTime)
        plt.text(0,0,' Time:\n ' + str(currentTime)  + ' s')
        f.savefig('Time_' + str(intSec) + "_" + str(intDecimal))
        plt.close(f)
    #End of big.matrix for loop.
    ##big.matrix <- big.matrix[1:move.count,]
    #### End Standard Execution

np.savetxt("vel1_x.csv", vel1_x, delimiter=",")
np.savetxt("vel2_x.csv", vel2_x, delimiter=",")
np.savetxt("pos1_x.csv", pos1_x, delimiter=",")
np.savetxt("pos2_x.csv", pos2_x, delimiter=",")
np.savetxt("std1.csv", std1, delimiter=",")
np.savetxt("std2.csv", std2, delimiter=",")
np.savetxt("resolution.csv", plot_res, delimiter=",")
np.savetxt("timeOutput.csv", time, delimiter=",")
#write.csv(peak_width,file="peak_width.csv",row.names=FALSE)
#write.csv(plot_res,file="plot_res.csv",row.names=FALSE)
##### End picture drawing
