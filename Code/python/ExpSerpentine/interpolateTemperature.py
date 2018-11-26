import numpy as np
import matplotlib.pyplot as plt

class gradientData:
    def __init__(self,x,t,T):
        self.x = x
        self.t = t
        self.T = T
        self.tstep = t[1] - t[0]
        self.molec_index = np.array([])

    def interpolate(self,x1,x2,y1,y2,x):
        #return x
        if(x1 == x2):
            return y1
        m = 1.0*(y1-y2)/(x1-x2)
        b = y1 - m*x1
        y = m*x+b
        return y

    def getTimeIndices(self,t):
        t1 = 1.0*t/self.tstep
        if(t1 == int(t1)):
            t2 = int(t1)
        else:
            t2 = int(t1 + 1)
        t1 = int(t1)
        return (t1,t2)

    def getPositionIndices(self,x):
        length = len(self.x)
        xmax = self.x[length-1]
        ind = int((x/xmax)*(length-1))
        ind2 = ind
        while(True):
            if(x > self.x[ind]):
                if(ind + 1 == len(self.x)):
                    ind2 = ind
                    break
                elif(x < self.x[ind+1]):
                    ind2 = ind + 1
                    break
                else:
                    ind = ind + 1
            elif(x < self.x[ind]):
                if(ind == 0):
                    ind2 = ind
                    break
                else:
                    ind = ind - 1
            else:
                ind2 = ind
                break
        return (ind,ind2)

    def getTemperature(self,x,t):
        (x1,x2) = self.getPositionIndices(x)
        (t1,t2) = self.getTimeIndices(t)
        T1 = self.interpolate(self.x[x1],self.x[x2],self.T[t1][x1],self.T[t1][x2],x)
        T2 = self.interpolate(self.x[x1],self.x[x2],self.T[t2][x1],self.T[t2][x2],x)
        T = self.interpolate(self.t[t1],self.t[t2],T1,T2,t)
        return T

    def plotTemperature(self,x,t):
        (x1,x2) = self.getPositionIndices(x)
        (t1,t2) = self.getTimeIndices(t)
        T1 = self.interpolate(self.x[x1],self.x[x2],self.T[t1][x1],self.T[t1][x2],x)
        T2 = self.interpolate(self.x[x1],self.x[x2],self.T[t2][x1],self.T[t2][x2],x)
        T = self.interpolate(self.t[t1],self.t[t2],T1,T2,t)
        plt.plot(self.x,self.T[t1],label = ("Time: ",  self.t[t1]))
        plt.plot(self.x,self.T[t2],label = ("Time: ", self.t[t2]))
        plt.xlabel("Position")
        plt.ylabel("Temperature")
        plt.scatter(x,T,label = ("x:",x," T:",round(T,2),"t:",t))
        plt.legend(loc = 'upper left')
        plt.title("Interpolation Temperature")
        plt.show()
        return T

# x_molecules = np.array([1.0,2.3,4.5,5.67,7.1])
# t = np.array([0,1,2,3,4,5,6,7,8,9,10])
# x = np.array([0,2,4,6,8,10,12,14,16,18,20])
# temp = np.array([10,10.2,10.8,11.5,11.7,12.2,13.1,13.3,13.5,13.7,15.0])
# T = np.array([temp])
# for i in range(0,11):
#     for j in range(len(temp)):
#         temp[j] = temp[j] + np.random.random()
#     T = np.append(T,[temp],axis=0)
# np.savetxt("T.csv", T, delimiter=",")
# T = np.loadtxt('temperature.csv',delimiter = ',')
# x = np.loadtxt('colPosition.csv',delimiter = ',')
# t = np.loadtxt('time.csv',delimiter = ',')
#
# data = gradientData(x,t,T)
# data.plotTemperature(1.2,2.5)
# data.plotTemperature(1.2,7.4)
# data.plotTemperature(0,2.5)
# data.plotTemperature(20,2.5)
# data.plotTemperature(1.2,0)
# data.plotTemperature(1.2,10)
# data.plotTemperature(13,6.2)
# data.plotTemperature(8.5,8.3)
