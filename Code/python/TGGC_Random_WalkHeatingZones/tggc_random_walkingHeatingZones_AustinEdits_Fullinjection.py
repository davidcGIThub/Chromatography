import numpy as np
import pandas as pd
import random as rnd
import matplotlib.pyplot as plt
import HelperFunctionsR as hfr
import math
import seaborn as sns
import time

# Defines colors and plot markers for each compound
# But my new plotting version doesn't even use it
colors = np.array(["red","green","magenta","black","blue","gray"])
markers = np.array(["o","^","P","x","D","v","s"])

j = 11.0 #number of compounds

n_mol = 50.0 # Number of molecules per compound

sig_sq = 4.073e-6 #3.792*exp(-4) #1*exp(-6)    # Involved in computing the diffusion term D_t
                                                # although I don't have a clear understanding
                                                # of the role it plays in defining D_t

taylor = 1.647*10**(-1)             # Same as above, it does something for D_t but I don't
                                    # know what

X = np.array([10,12,14,16,18,20,22,24,26,28,30]) #alkane length
C =  - 0.47406*X - 7.22173          # used in computing k (retention factor) but I don't know
                                    # specifically what it represents
h = 2284 + 819*X                    # Same as C. It does something in computing k, but I
                                    # don't know what

delta_t = 0.0001 # time step: number of n's per second
R = 1.987 # Boltzman's constant

u_m = 55.0 # velocity of mobile phase [cm/sec]

run_time = 1500000000

#---------------------------------------------------------------

# Read in simulated temperature data
# > Note to David: You'll have to change the filePath
# >                 to wherever you save the data table files
# >                 on your machine. The files are on the J
# >                 drive in the same folder as the other
# >                 simulated gradients. Just look for
# >                 "time.csv", "postition.csv", and
# >                 "temperature.csv"

t1 = time.time()
filePath = r'J:\groups\fluxlab\TGGC\Austin\Code Developed\2018_8_14 Tolley TGGC Model Python'
# Note: If you access these files directly from the J drive it takes
#       7 times as long to read in the data as if you have it saved
#       locally on your machine. Right now it takes 18 seconds on my
#       computer

timeArr = np.genfromtxt(filePath+r'\time.csv',delimiter=',')
posArr = np.genfromtxt(filePath+r'\position.csv',delimiter=',')
tempArr = np.genfromtxt(filePath+r'\temperature.csv',delimiter=',')

##timeArr = np.genfromtxt(filePath+r'\timeSmooth.csv',delimiter=',')
##posArr = np.genfromtxt(filePath+r'\positionSmooth.csv',delimiter=',')
##tempArr = np.genfromtxt(filePath+r'\temperatureSmooth.csv',delimiter=',')

t2 = time.time()
print('File Read Time: ' + str(t2 - t1))

col_length = max(posArr)*100

#---------------------------------------------------------------

##### Start Standard Execution
detector = np.zeros((int(j),int(n_mol)))
x = np.tile(np.linspace(1,n_mol,int(n_mol))/10 , (int(j),1) )
np.apply_along_axis(np.random.shuffle,-1,x)
t2 = time.time()
print("x creation: " , t2 - t1)

#---------------------------------------------------------------

# Vectorized temperature reading test
t1 = time.time()
n = 0
timeIn = 108.9
delta_t_Sim = 1

if timeIn < 1:

    Tout_1 = np.interp(x,posArr,tempArr[0,:])

else:

    ind1 = np.argmin(abs(timeArr - timeIn))
    if timeArr[ind1] > timeIn:
        ind = [ind1-1,ind1]
    else:
        ind = [ind1,ind1+1]
    delta_Temp = tempArr[ind[1],:] - tempArr[ind[0],:]
    tempInterpIntermediate = delta_Temp*(timeIn - int(timeIn))/delta_t_Sim + tempArr[int(timeIn)-1,:]
    TempOut = np.interp(x,posArr,tempInterpIntermediate)

t2 = time.time()
print('Temp Inpterp Time: ',t2 - t1)

#---------------------------------------------------------------

# Extreme test gradient
T_gradTest = np.linspace(500,100,len(posArr))

#---------------------------------------------------------------

res_end = np.zeros((int(j),2))
eluted_1 = 1
eluted_2 = 1
detected_1 = np.array([])
detected_2 = 0
detection_AF = []   # "AF" stands for Austin Foster
                    # I don't want you getting the
                    # wrong idea....

t = time.time()
tStat = []
meanStat1 = []
stdStat1 = []
meanStat2 = []
stdStat2 = []
for n in range(0,run_time):

    # Austin's temperature interpolation method
    timeIn = n*delta_t
    delta_t_Sim = 1

    if timeIn < 1:

        tempInterpIntermediate = tempArr[0,:]
        TempOut = np.interp(x,posArr*100,tempArr[0,:])

    elif timeIn >= timeArr[-1]:

        tempInterpIntermediate = tempArr[-1,:]
        TempOut = np.interp(x,posArr*100,tempArr[-1,:])

    else:

        ind1 = np.argmin(abs(timeArr - timeIn))
        if timeArr[ind1] > timeIn:
            ind = [ind1-1,ind1]
        else:
            ind = [ind1,ind1+1]
        delta_Temp = tempArr[ind[1],:] - tempArr[ind[0],:]
        tempInterpIntermediate = delta_Temp*(timeIn - timeArr[ind[0]])/delta_t_Sim + tempArr[ind[0],:]
        TempOut = np.interp(x,posArr*100,tempInterpIntermediate)

##    tempInterpIntermediate = tempArr[0,:]
##    TempOut = np.interp(x,posArr*100,tempArr[0,:])

##    TempOut = np.interp(x,posArr*100,T_gradTest)

    T_all = TempOut.flatten('F')

    
    # Dr. Tolley's transport computation method.
    # I didn't change anything here
    C_all = C
    tmp_h = np.repeat([h],len(T_all)/len(h),axis = 0).flatten()
    tmp_C_all = np.repeat([C_all],len(T_all)/len(C_all),axis=0).flatten()
    k_all = np.exp(tmp_h/(R*T_all) + tmp_C_all)
    D_t = ( (T_all**1.75)*sig_sq + ((1+6*k_all+11*k_all**2)/((1+k_all)**2))*(u_m**2)*taylor/(T_all**1.75) )*delta_t
    x = np.reshape(x.flatten('F') + (u_m*delta_t + hfr.pyRNorm(j*n_mol,0,(2*D_t)**(.5))) * ((np.random.rand(int(j*n_mol)) < (1.0/(1+k_all))) + 0) , np.shape(x),'F') #adding the zero makes it numeric


    x[x<0] = -100       # This is part of the workaround I used to be able to detect
                        # molecules and then make sure they wouldn't be detected
                        # twice without having to do weird things to the x array.
                        # All I do is when a molecule exits the column I say that
                        # it has a negative position and this line makes sure that
                        # it always keeps a negative position

##    # Compute mean and standard deviation
##    if(n%1000 == 0):
##        tStat.append(n*delta_t)
##        xComp1 = x[0,x[0,:]>=0]
##        xComp2 = x[1,x[1,:]>=0]
##        if len(xComp1)==0:
##            stdStat1.append(stdStat1[-1])
##            meanStat1.append(meanStat1[-1])
##        else:
##            stdStat1.append(np.std(xComp1))
##            meanStat1.append(np.mean(xComp1))
##
##        if len(xComp2)==0:
##            stdStat2.append(stdStat2[-1])
##            meanStat2.append(meanStat2[-1])
##        else:
##            stdStat2.append(np.std(xComp2))
##            meanStat2.append(np.mean(xComp2))

    # New method for checking for detection
    if (max(x.flatten('F'))>col_length):
        for k in range(int(j)):
            all_detected = np.where(x[k,:] > col_length)[0]
            add = [n*delta_t]*len(all_detected)
            detection_AF.extend(add)
            x[k,all_detected] = -100
##        all_detected_1 = np.where(x[0,:] > col_length)[0]
##        all_detected_2 = np.where(x[1,:] > col_length)[0]
##        add_1 = [n*delta_t]*len(all_detected_1)
##        add_2 = [n*delta_t]*len(all_detected_2)
##        detection_AF.extend(add_1)
##        detection_AF.extend(add_2)
##        x[0,all_detected_1] = -100
##        x[1,all_detected_2] = -100
    if (max(x.flatten('F'))<0):
        break

##    # Plotter            
##    if(n%1000000 == 0):
##        print("pass: " , t - time.time())
##        t = time.time()
##
##        # Austin Plot
##        f,(ax1,ax2) = plt.subplots(2, sharex=True)
##        ax1.hist(x.flatten('F'))
####        ax1.hist(x[1,:])
##        ax2.plot(posArr*100,tempInterpIntermediate)
####        ax2.plot(posArr*100,T_gradTest)
##        ax1.set_title('Time: ' + str(n*delta_t) + ' s')
##        plt.show()

##        # mean and std Plot
##        f1,(ax3,ax4) = plt.subplots(2, sharex=True)
##        ax3.plot(tStat,meanStat1,label='Peak 1')
##        ax3.plot(tStat,meanStat2,label='Peak 2')
##        ax4.plot(tStat,stdStat1)
##        ax4.plot(tStat,stdStat2)
##        ax3.set_title('Time: ' + str(n*delta_t) + ' s')
##        ax3.set_ylabel('Mean')
##        ax4.set_ylabel('Standard Deviation')
##        ax3.legend(loc=0)
##        ax3.grid(b=1)
##        ax4.grid(b=1)
##        plt.show()
##
##        # T and k distribution plots
##        ax3.clear()
##        ax4.clear()
##        f1,(ax3,ax4) = plt.subplots(2, sharex=False)
##        ax3.hist(T_all)
##        ax4.hist(k_all)
##        ax3.set_title('Time: ' + str(n*delta_t) + ' s')
##        ax3.set_ylabel('T Dist')
##        ax4.set_ylabel('k Dist')
##        plt.show()

        
t3 = time.time()
print("loop: " , t3 - t2)

# Detection Plotter
plt.hist(np.asarray(detection_AF)/60,bins = 300,range=(0,max(detection_AF)/60+2))
plt.xlabel('Elution Time (min)')
plt.show()

np.savetxt('detectionData.csv',detection_AF,delimiter=',')

### mean and std Plot
##f1,(ax3,ax4) = plt.subplots(2, sharex=True)
##ax3.plot(tStat,meanStat1,label='Peak 1')
##ax3.plot(tStat,meanStat2,label='Peak 2')
##ax4.plot(tStat,stdStat1)
##ax4.plot(tStat,stdStat2)
##ax3.set_title('Time: ' + str(n*delta_t) + ' s')
##ax3.set_ylabel('Mean')
##ax4.set_ylabel('Standard Deviation')
##ax3.legend(loc=0)
##ax3.grid(b=1)
##ax4.grid(b=1)
##plt.show()
