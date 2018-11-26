import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns

########################
###   NEW EQUATION:  This is a step by step program to generate a random walk model
###     for a varietly of moving thermal gradient profiles.
###		This code was UPDATED in Riva with Samuel in May 2016
####################
j = 2 # number of compounds
nMol = 100 #****Number of molecules per compound used in simulation
colors  = ['red','green','blue','grey30']
markers = np.array(["o","^","P","x","D","v","s"])
sigmaNot = 10000*math.sqrt(1.0/12)
gamma1 = 5.113e-3 # This is the molecular diffusion coefficient when using both temp and pressure.
taylor = 1.647e-1  #  This is the most recent based on solving Anzi's data, using R=1.987.  (141105)
                   # This is the dispersion coefficient due to resistence to transport and taylor dispersion.
gamma3 = 7.676e-8 # This is resistence to flow for adsorption/deporption (no pressure variable)
gamma2 = 1.217e-9 # This is Golay's formula for resistence to transport when using both temp and pressure in diffusion expression.
chkSm = 0
#C <- c(-12.5075, -12.7202)
#h <- c(11416, 11521)
X = np.array([9.0,10.0]) #alkane length
C = -7.22173 - 0.47406*X  #This is the formula for calculating C where X is the number of carbons.  Includes phase ratio
h = 2284.0 + 819*X    #  calculation of enthalpy (need to check this)
h = h* 4.18442979  # This is to change enthalpy from kcal units to Joules in the units.
C = np.array([-11.1277,-12.4829]) #???? Why change the values?
h = np.array([34926.27,40845.71])
C = np.array([-11.1277,-11.2829]) #???? Why change the values?
h = np.array([34926.27,35245.71])
delta_t = 0.0001 # 0.001 ***# time step: Reciprocol is the number of time steps per second
R = 8.3144621 #  Boltzman's constant in Joules per Kelvin mole.  I used to use R = 1.987 kcal per Kelvin mole.
colLength = 20.0#0 #*** column length [cm]
colLength = 1000.0 # measured in cm
T0 = 323.0 # initial temperature [Kelvin, K]
T0 = 313.0 #???? Why change the values?
T0 = 290.0
diameter = 0.0100 #inner diameter of the column in centimeters.
#diameter <- 0.0001 #measured in meters for test.
p_i = 57.16  #This is the column inlet pressure by gage measured in psi.  For true inlet must add atmosphere pressure.
#p.i <- 30
p_o = 14.69595   #This is the outlet pressure or Atmospheric pressure in psi.
p_o = 13.8 #???? Why change the values?
p_i = p_i*6894.757 # This converts pressure from psi units to Pascals  (Anzi solution)
p_o = p_o*6894.757 # This converts pressure from psi units to Pascals
#p.i <- p.i*10000/1.450377  #converting psi to Pa (pascals)  (Wiki Solution)
#p.o <- p.o*10000/1.450377  #converting psi to Pa
p_i = p_i + p_o  #This is the true inlet pressure.
w = 6.0 # velocity of gradient [cm/sec]
tpRate = 0.5 #degrees C per second. This is the rate at which the column is heated after a gradient is established.
#tp.rate <- 1.7
tpRate = .25 #This is the more reasonable rate of raising temperature under TPGC
#             rep(0,times=1e4*(col.length-ramp.length) )
TempGranularity = 1.0e2#4 # granularity of temperatures across column. data points/cm
powerAce = 100.0


###################
#
#		Here is the temperature profile function.
#		It produces the initial heat gradient profile that we heat.
#		We need to make this a function call to either a data file or a function.
#
#		The values needed are:
#			Temp.granularity
#			ramp.length
#			ramp.height
#			zone.len
#			window
#			col.length
#
#		Return a vector  Temp.vec that has the profile of temperature along
#		the column.  We will then heat this profile up.
#
##################
plateLen =  10.0 #60 #cm  Note that the plate width is calculated by subtracting the plate.len
                  #from the column length and then dividing the difference by the number of
                  # turns.  This are called isotherms.
                  #
decayConst = 0.5 #0.05  # THis is the exponential decay factor, giving the rate of decay
TRamp = 50.0 #height of the temperature ramp; delta T across plate.  This is not T.0
expoDecay = np.exp(-np.arange(1,(plateLen*TempGranularity)+1)/TempGranularity*decayConst)*TRamp
     #This is determines the exponential decay
     #
#plot(expo.decay[(1:plate.len)*Temp.granularity],type='l', xlab='Temperature', main=paste("PlateLength: ",plate.len, " DecayConst: ", decay.const, "NumbCorners: ", n.corners))
# Execute the above line to look at exponential temperature decay down plate.
nCorners = 200# number of switchbacks, and rungs  #This turns out to be the number of isotherms.
# n.corners <- col.length-plate.len
lenCross = (colLength-plateLen)/nCorners  #This gives the length of each isotherm.
seqStart = 1
cornerSpacing = 1.0*(len(expoDecay)-seqStart)/(nCorners) ###check
corner = np.round(np.arange(seqStart,len(expoDecay)+1, cornerSpacing)) #indexes the position of each of the corners
            # on the exponential decay curve.
TempVec = np.array([])
for jj in range(0,nCorners):
    #This gives the temperature at each of the corners
    TempVec = np.append(TempVec,expoDecay[int(corner[jj]-1):int(corner[jj+1])])
    #  This determines the temperature at each point along the isotherm.
    TempVec = np.append(TempVec, np.repeat(expoDecay[int(corner[jj+1]-1)],lenCross*TempGranularity))
TempVec = TempVec[0:int(colLength*TempGranularity)]
#plot(Temp.vec[1:(col.length*10)*Temp.granularity/10],type='l', xlab='Temperature', main=expression(paste("PlateLength: ",plate.len, " DecayConst: ", decay.const, "NumbCorners: ", n.corners))) # Execute to look at temperature along column
window = 100   #This is for the ksmooth function, below.  It should be proportional to the
                # temperature granularity above.  Change either and you must change the other.
#Temp.vec <- ksmooth(1:length(Temp.vec),Temp.vec,'normal',window)$y  REMOVED #This smooths the abrupt step function produced above.
#  NOTE!!!  ksmooth takes a long time.  Samuel is not sure it is necessary.   I am thinking that if
#  it turns out to be necessary, then we could archive these gradient curves and recall as necessary.
#
###Temp.vec <- unlist(Temp.vec) REMOVED  #Here we make the 'list' structure from the object above into an array.
#plot(Temp.vec[(1:(length(Temp.vec)/1e3))*1e3],type='l')
#******************
#  Here I set the gradient to zero, giving programmed temperature or isothermal
#
#******************
#Temp.vec <-rep(0,length(Temp.vec))
#tp.rate <- 0
###################
#
#		This is the end of the temperature profile function
#
###################
#address the temperature /100, get the temp raised to the power
Temp2PowVec = np.arange(0,1000+1/powerAce,1/powerAce)**1.646
#address the temperature /100, get the temp. No Power here (HDT)
TempNoPowVec = np.arange(0,1000+1/powerAce,1/powerAce)
#################
#
#	The vector 'TempNoPow.vec' containing the basic temperature profile or template
#		is returned from temperature profile function call.
#	The vector 'Temp2Pow' is sequence of temperatures raised to the 1.646 power
#		to adjust for the viscosity of helium.
#
#################


runTime = 1500000  #Number of potential steps to get compounds out of column.
#run.time <- 6000
pauseCount = 10000  #number of frames to skip in plotting results.  The 'big.matrix'
                   # takes a snapshot at every 'pause.count' time for illustrating the
                   # progress of the molecules.
bigMatrix = np.zeros((runTime/pauseCount+1,nMol*j))
#  This has the location of every molecule at each recorded time point.  Only
#   time points recorded are every "pause.count" apart.
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
moveCount = 0
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

Tmax = delta_t * runTime * tpRate + TRamp + T0
Tmin = T0
Tdelta = Tmax-Tmin
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
    TempPlus = TempVec + tpRate*n*delta_t + T0
    #print(TempPlus)
    #print(len(TempPlus))
       #Here we heat up the whole column an amount 'tp.rate*delta.t
           # from the last time increment.  To avoid numerical problems
           # we simply multiply the temperature increment amount by the
           # number of time steps, n, and add to initial temperature.
           # Initial temperature is T.0 plus the temperature gradient profile.
       #************
       # Now I hold the temperature constant.  This is for isothermal of gradient
       #  NOTE: Comment this out if allow heating of column over time.
       #Temp.plus <- Temp.vec + T.0
       #
       #*************
    temp = abs(x*TempGranularity+1)
    temp[temp > colLength*TempGranularity] = colLength*TempGranularity
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
    velocity_x = ((p_i**2-p_o**2)*TempPlus[temp.astype(int)] * (diameter/2)**2 / ((16*HeFactor*totInt )*np.sqrt(abs(p_i**2-(p_i**2-p_o**2)*TSumPow[temp.astype(int)]/totInt ))))
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
    # mu1 = mean(x[0,])
    # mu2 = mean(x[1,])
    # if mu1 <= colLength:
    #     totalTime1 = n*delta_t
    #     speed1 = velocity_x[0]/(1+k_all[0])
    #     spread1 = np.std(x[0,])*(1+k.all[0])/(60*velocity_x[0])
    # if mu2 <= colLength:
    #     totalTime2 = n*delta_t
    #     speed2 = velocity_x[1]/(1+k_all[1])
    #     spread2 = np.std(x[1,])*(1+k_all[1])/(60*velocity_x[1])
    # if mu1 <= 20:
    #    time1_20 <- n*delta.t
    #    speed.1.20 <- velocity.x[1]/(1+k.all[1])
    #    spread.1.20 <- sd(x[1,])*(1+k.all[1])/(60*velocity.x[1])
    #  }
    #  if(mu.2<= 20)
    #  {
    #    time.2.20 <- n*delta.t
    #    speed.2.20 <- velocity.x[2]/(1+k.all[2])
    #    spread.2.20 <- sd(x[2,])*(1+k.all[2])/(60*velocity.x[2])
    #  }
    #  if(mu.1<= col.length/4)
    #  {
    #    time.1.fourth <- n*delta.t
    #    speed.1.fourth <- velocity.x[1]/(1+k.all[1])
    #    spread.1.fourth <- sd(x[1,])*(1+k.all[1])/(60*velocity.x[1])
    #  }
    #  if(mu.2<= col.length/4)
    #  {
    #    time.2.fourth <- n*delta.t
    #    speed.2.fourth <- velocity.x[2]/(1+k.all[2])
    #    spread.2.fourth <- sd(x[2,])*(1+k.all[2])/(60*velocity.x[2])
    #  }
    #  if(mu.1<= col.length/2)
    #  {
    #    time.1.half <- n*delta.t
    #    speed.1.half <- velocity.x[1]/(1+k.all[1])
    #    spread.1.half <- sd(x[1,])*(1+k.all[1])/(60*velocity.x[1])
    #  }
    #  if(mu.2<= col.length/2)
    #  {
    #    time.2.half <- n*delta.t
    #    speed.2.half <- velocity.x[2]/(1+k.all[2])
    #    spread.2.half <- sd(x[2,])*(1+k.all[2])/(60*velocity.x[2])
    #  }
     #  1/(1+exp( h/(R*( T.all )) + C.all )
    if n%500==0:
        if np.max(x[0,]) > colLength:
            allDetected1 = x[0,][x[0,] > colLength]
            newDetected1 = allDetected1[allDetected1 != detected1]
            detected1 = np.concatenate(detected1, newDetected1)
            detector[0, newDetected1] = n*delta_t
        if np.max(x[1,]) > colLength:
            allDetected2 = x[1,][x[1,] > colLength]
            newDetected2 = allDetected2[allDetected2 != detected2]
            detected2 = np.concatenate(detected2, newDetected2)
            detector[1, newDetected2] = n*delta_t
            if len(detected_2)+len(detected_1)-2 >= 2*nMol:
                #break
                break
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
        bigMatrix[moveCount,] = x.flatten()
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
        #x = bigMatrix[moveCount,].reshape(j,nMol)
        for i in range(0,j):
           lab = round(4*np.std(x[i,],axis=0,ddof=1) , 2)
           ax1.scatter(x[i,], 1.0*np.arange(1,nMol+1)/nMol*Tdelta+Tmin,color = colors[i],s = 0.5,  marker = markers[i],label =lab)
        ax1.legend()
        tmp = (np.std(x[0,],axis=0,ddof=1)+np.std(x[1,],axis=0,ddof=1))/2.0
        peak_width = np.append(peak_width, 4*tmp)
        res = (np.mean(x[0,])-np.mean(x[1,]))/(4.0*tmp)
        plot_res = np.append(plot_res, res)
        xx = np.linspace(0,colLength-1,(colLength-1)/.1 + 1)
        yy = (TempVec[(xx*TempGranularity).astype(int)]+tpRate*delta_t*(moveCount)*pauseCount+T0)
        ax1.plot(xx,yy,color = 'red')
        #abline(h=( (delta.t*mm*pause.count*b) )/200,col='red') #programmed temperature
        ax2.set_ylabel('Density')
        for i in range(0,j):
           sns.kdeplot(x[i,],color=colors[i],ax=ax2) #density plot
        ##added DC #######################################
        #currentTime = pauseCount*delta_t*moveCount
        time = np.append(time,moveCount)
        plt.text(0,0,' Time:\n ' + str(moveCount)  + ' s')
        f.savefig('Time:' + str(moveCount))
        plt.close(f)
        moveCount = moveCount + 1
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
np.savetxt("time.csv", time, delimiter=",")
#write.csv(peak_width,file="peak_width.csv",row.names=FALSE)
#write.csv(plot_res,file="plot_res.csv",row.names=FALSE)
##### End picture drawing



# retention.time <- rowMeans(detector)
# stand.dev <- apply(detector,1,sd)
# resolutions <- rep(0,j)
#      for(i in 2:j){resolutions[i] <- (retention.time[i]-retention.time[i-1])/(2*(stand.dev[i]+stand.dev[i-1]))}
# (rbind(stand.dev,retention.time, resolutions))
# ##### Simulation Results
# plate.len <-  30 #cm
# decay.const <- 0.05
# T.ramp <- 90 #height of the temperature ramp; delta T across plate
# n.corners <- 20
# #stand.dev       0.02365578   0.02487849
# #retention.time 40.46400000  51.47850000
# # resolutions     0.00000000 113.47136009
# plate.len <-  30 #cm
# decay.const <- 0.1
# # stand.dev       0.02480225   0.01305582
# # retention.time 52.17100000  70.30250000
# # resolutions     0.00000000 239.46676175
# ## Double the length, but same number of turns shold make the gradient portions longer and isothermal shorter
#  plate.len <-  60 #cm
# decay.const <- 0.1  #first half same as plate.len=30, decay.const=0.1, but extends to 60 cm
# # stand.dev       0.0223098   0.01507557
# # retention.time 64.5135000  88.05500000
# # resolutions     0.0000000 314.84910248
# plate.len <-  60 #cm
# decay.const <- 0.05  # Same shape as plate.len= 30, decay.const= 0.1
# #stand.dev       0.02480225 0.008689364
# #retention.time 52.12900000 70.24950
# #resolutions     0.00000000 270.5230
# plate.len <- 120 #cm
# decay.const <- 0.025 # same shape as 30;0.1, and 60,0.05.
# #stand.dev       0.01882938   0.01930615
# # retention.time 52.05700000  70.14100000
# # resolutions     0.00000000 237.10176732
