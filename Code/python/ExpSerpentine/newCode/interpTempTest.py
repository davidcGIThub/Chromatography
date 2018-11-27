import interpolateTemperature as interpT
import numpy as np

x_molecules = np.array([1.0,2.3,4.5,5.67,7.1])
t = np.array([0,1,2,3,4,5,6,7,8,9,10])
x = np.array([0,2,4,6,8,10,12,14,16,18,20])
temp = np.array([10,10.2,10.8,11.5,11.7,12.2,13.1,13.3,13.5,13.7,15.0])
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
print(data.getGradient(2.1))
data.plotTemperature(1.2,2.5)
data.plotTemperature(1.2,7.4)
data.plotTemperature(0,2.5)
data.plotTemperature(20,2.5)
data.plotTemperature(1.2,0)
data.plotTemperature(1.2,10)
data.plotTemperature(13,6.2)
data.plotTemperature(8.5,8.3)
