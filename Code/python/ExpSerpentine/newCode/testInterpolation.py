import numpy as np
import interpolateTemperature2 as interpT


t = np.array([1,2,3,4,5,6,7,8,9,10])
x = np.array([2,4,6,8,10,12,14,16,18,20])
temp = np.array([10.2,10.8,11.5,11.7,12.2,13.1,13.3,13.5,13.7,15.0])
T = np.array([temp])
for i in range(0,11):
    for j in range(len(temp)):
        temp[j] = temp[j] + np.random.random()
    T = np.append(T,[temp],axis=0)
#np.savetxt("T.csv", T, delimiter=",")
#T = np.loadtxt('temperature.csv',delimiter = ',')
#x = np.loadtxt('colPosition.csv',delimiter = ',')
#t = np.loadtxt('time.csv',delimiter = ',')

data = interpT.GradientData(x,t,T)
data.plotTemperature(6.7,100)
data.plotTemperature(0,5.0)
data.plotTemperature(9,5.0)
data.plotTemperature(3.6,5.0)
data.plotTemperature(8.5,5.0)
