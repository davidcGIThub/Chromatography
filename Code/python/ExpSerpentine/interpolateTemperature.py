import numpy as np
import matplotlib.pyplot as plt


def interpolateTemp(x1,x2,T1,T2,x):
    #return x
    m = 1.0*(T1-T2)/(x1-x2)
    b = T1 - m*x1
    T = m*x+b
    return T

#tData 2X2 array col 1 = pos data, col 2 = Temp data, x = pos of desired data
def interpolateTempTime(t1,t2,t1Data,t2Data,x,t):
    T1 = interpolateTemp(t1Data[0][0],t1Data[0][1],t1Data[1][0],t1Data[1][1],x)
    T2 = interpolateTemp(t2Data[0][0],t2Data[0][1],t2Data[1][0],t2Data[1][1],x)
    T = interpolateTemp(t1,t2,T1,T2,t)
    return T

t = np.array([0,1,2,3,4,5,6,7,8,9,10])
x = np.array([0,2,4,6,8,10,12,14,16,18,20])
temp = np.array([10,10.2,10.8,11.5,11.7,12.2,13.1,13.3,13.5,13.7,15.0])
T = np.array([temp])
for i in range(0,11):
    for j in range(len(temp)):
        temp[j] = temp[j] + np.random.random()
    T = np.append(T,[temp],axis=0)
print(T)

time = 1.3
x = 3

plt.plot(x,T[0],x,T[1],)
plt.show()
