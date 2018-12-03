import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import Parameters_ExpSerpentine as p
import interpolateTemperature2 as interpT

########################
###   NEW EQUATION:  This is a step by step program to generate a random walk model
###     for a varietly of moving thermal gradient profiles.
###		This code was UPDATED in Riva with Samuel in May 2016
####################
TemperatureGradients = np.loadtxt(p.gradientFileName,delimiter = ',')
colPosition = np.loadtxt(p.positionFileName,delimiter = ',') * p.posUnits
timeData = np.loadtxt(p.timeFileName,delimiter = ',')
GData = interpT.GradientData(colPosition,timeData,TemperatureGradients)
colLength = colPosition[len(colPosition)-1] - colPosition[0]  # Column length [cm]

#C = -7.22173 - 0.47406*p.X                # Formula for calculating C where X is the number of carbons.  Includes phase ratio
C = -9.7309 - 0.1552*p.X             #Formula using hard coded C values
#h = (2284.0 + 819*p.X) * 4.18442979       #  Calculation of enthalpy (Joules)
h = (7659.66012301 + 76.4301505179*p.X) * 4.18442979       #  Calculation of enthalpy (Joules) using hard code values to calc

p_i = p.p_i*6894.757          # This converts pressure from psi units to Pascals
p_o = p.p_o*6894.757          # This converts pressure from psi units to Pascals
p_i = p_i + p_o             #This is the true inlet pressure.

Temp2PowVec = np.arange(0,1000+1/p.powerAce,1/p.powerAce)**1.646        #address the temperature /100, get the temp raised to the power
#################
#	The vector 'Temp2Pow' is sequence of temperatures raised to the 1.646 power
#		to adjust for the viscosity of helium.
#################
runTime = int( (timeData[len(timeData)-1] - timeData[0]) / p.delta_t)  #Number of potential steps to get compounds out of column.
pauseCount = int(1/p.delta_t) #number of frames to skip in plotting results.  The 'big.matrix'
                            # takes a snapshot at every 'pause.count' time for illustrating the
                            # progress of the molecules.
##### Start Standard Execution
x = np.zeros((p.j,p.nMol))
                                    # This matrix keeps track of x-axis location of each
                                    # molecule of each analyte. Initially each is at point
                                    # '0' as indicated here.  The y-axis location is simply the
                                    # scaled value of the molecule's index number in the matrix x.
detector = np.zeros((p.j,p.nMol))
moleculeExitTime = np.zeros((p.j,p.nMol))

for i in range(0,p.j):
  x[i,] = np.random.choice(p.nMol, p.nMol, replace=False)*1.0*p.injectionWidth/p.nMol
# The above loop returns an order from 1 to n.mol for each molecule.  Since it is divided
# by n.mol, in essence we get each molecule uniquely with a value between 0 and 1.
# Multiplying by injection.width gives a plug that has the n.mol molecules uniformly distributed
# from 0 cm to injection.width cm at the beginning of the column.
#   Remember, x[i,] is the vector of locations in x-axis of i-th analyte.  This location is used
#	at each step for each molecule to determine temperature and velocity.
#	This is also used to plot location of the molecules.  The y-axis values for plot are determined
#	from the matrix column index.  (Rows indicate analyte and column indicate individual molecules.)
currentTime = timeData[0]
plot_res = np.zeros((p.j-1,0))
#added_dc#############################################################
vel_x = np.zeros((p.j,0))
pos_x = np.zeros((p.j,0))
std = np.zeros((p.j,0))
time = np.array([])
Tmax = np.amax(TemperatureGradients)
Tmin = np.amin(TemperatureGradients)
Tdelta = Tmax-Tmin

for n in range(0,runTime+1):
#for n in range(0,2):
    currentTime = round(n*p.delta_t,3)
    TempPlus = GData.getTemperatures(n*p.delta_t,p.tempGranularity)
    temp = abs(x*p.tempGranularity+1)
    temp[temp >= len(TempPlus)] = len(TempPlus) - 1
    T_all = TempPlus[temp.astype(int)].flatten()
    C_all = C
    tmp_C_all = np.repeat([C_all],len(T_all)/len(C_all),axis=0).flatten('F')
    tmp_h = np.repeat([h],len(T_all)/len(h),axis = 0).flatten('F')
    k_all = np.exp(tmp_h/(p.R*T_all)+tmp_C_all)
       # Here I adjust the diffusion calculation according to the position of each molecule and the temperature
       #   at that location, taking adjusting the velocity resulting from the different pressure and viscosity values.  (HDT)
       # 	To do this there are several steps.
       #	First step is to determine a numerical value for the integral of temperature times viscosity.
    TSumPow = np.cumsum(Temp2PowVec[(TempPlus*p.powerAce).astype(int)-1]) / p.tempGranularity
       # integral of Temp^1.646
         # Note baseline viscosity cancels out.
         # T.sum.pow[length(T.sum.pow)]
         # is sum or integral from 0 to end.
    totInt = TSumPow[int(colLength*p.tempGranularity)-1]  #This is the integral from 0 to L of T^1.646.
    HeFactor = 0.000018662/(273.15**.646)  #Mult factor viscosity of He solvent.  This is delta.not in notes.
    CStar = (p_i**2-p_o**2)/(2*HeFactor*TSumPow[int(colLength*p.tempGranularity)-1])
         ## C.star is =to (pi^2-po^2)/(2*integral over L of T^1.646)
         #	I am doing the integral as simply the sum of the histogram pieces.  I could
         #	improve this step by using Simpson's rule or the Newton binomial trick if the
         #	Temperature can be modeled as a simple function of two variates. (HDT)
    velocity_x = ( (p_i**2-p_o**2)*TempPlus[temp.astype(int)] * (p.diameter/2)**2 / ((16*HeFactor*totInt )*np.sqrt(abs(p_i**2-(p_i**2-p_o**2)*TSumPow[temp.astype(int)]/totInt ))) )
    p_x = np.sqrt(abs((p_i**2-(p_i**2-p_o**2)*TSumPow[temp.astype(int)]/totInt)) ).flatten()
       #This is the pressure at x
    TempPrssRat = (T_all**1.75)/p_x
    Dt = 2*(( TempPrssRat*p.gamma1 + ((1+6*k_all+11*k_all**2)/((1+k_all)**2))*velocity_x.flatten()**2*p.gamma2/(TempPrssRat) + p.gamma3*velocity_x.flatten()**2*k_all/(1+k_all) )/(1+k_all))
    sigmaDelta = np.sqrt(abs(Dt*p.delta_t))
    x = np.reshape((x.flatten() + ( velocity_x.flatten()*p.delta_t/(1 + k_all) + np.random.normal(0,sigmaDelta,p.j*p.nMol))),np.shape(x))
    prevDetector = np.copy(detector)
    detector[np.where(x > colLength)] = 1
    moleculeExitTime[np.where(detector != prevDetector)] = currentTime
    if n % pauseCount == 0:
        if np.min(x) > colLength:
            break
        print(currentTime)
        time = np.append(time,currentTime)
        vel_x = np.append(vel_x,np.array([np.mean(velocity_x,axis = 1)]).T,axis=1)
        currentMeanPos = np.array([np.mean(x,axis = 1)]).T
        pos_x = np.append(pos_x,currentMeanPos,axis=1)
        std = np.append(std,np.array([np.std(x,axis=1,ddof=1)]).T,axis=1)
        f,(ax1,ax2) = plt.subplots(2, sharex=True)
        plt.xlim(0,colLength)
        plt.ylim(Tmin-Tmax)
        ax1.yaxis.set_label_position("left")
        ax2.yaxis.set_label_position("left")
        ax1.set_ylabel("Temperature (K)")
        for i in range(0,p.j):
           lab = str(i) + ": " + str(round(4*np.std(x[i,],axis=0,ddof=1) , 2))
           ax1.scatter(x[i,], 1.0*np.arange(1,p.nMol+1)/p.nMol*Tdelta+Tmin,color = p.colors[i],s = 0.5,  marker = p.markers[i],label =lab)
        ax1.legend(loc = 'best')
        tmp = np.array([np.std(x,axis=1,ddof=1)]).T
        tmp = (tmp[0:len(tmp)-1] + tmp[1:len(tmp)])/2
        res = (currentMeanPos[0:len(currentMeanPos)-1]-currentMeanPos[1:len(currentMeanPos)]) / (4.0*tmp)
        if(len(res) < 2):
            res = res.T
        plot_res = np.append(plot_res, res,axis=1)
        xx = np.linspace(0,colLength-1,(colLength-1)/.1 + 1)
        yy = (TempPlus[(xx*p.tempGranularity).astype(int)])
        ax1.plot(xx,yy,color = 'red')
        ax2.set_ylabel('Density')
        for i in range(0,p.j):
           sns.kdeplot(x[i,],color=p.colors[i],ax=ax2) #density plot
        intSec = int(currentTime)
        intDecimal = int((currentTime - intSec)*100)
        plt.text(0,0,' Time:\n ' + str(currentTime)  + ' s')
        f.savefig('Time_' + str(intSec) + "_" + str(intDecimal))
        plt.close(f)
    #### End Standard Execution
np.savetxt("molecExitTime.csv", moleculeExitTime, delimiter=",")
np.savetxt("vel_x.csv", vel_x, delimiter=",")
np.savetxt("pos_x.csv", pos_x, delimiter=",")
np.savetxt("std.csv", std, delimiter=",")
np.savetxt("resolution.csv", plot_res, delimiter=",")
np.savetxt("timeOutput.csv", time, delimiter=",")
##### End picture drawing