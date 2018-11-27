import numpy as np
import matplotlib.pyplot as plt

class GradientData:
    def __init__(self,x,t,T):
        self.position = x
        self.time = t
        self.Temperatures = T
        self.tstep = t[1] - t[0]

    def interpolateTempTime(self,t):
        (t1,t2) = self.getTimeIndices(t)
        T1 = self.Temperatures[t1]
        T2 = self.Temperatures[t2]
        if(t1 == t2):
            return T1
        m = 1.0*(T1-T2)/(self.time[t1]-self.time[t2])
        b = T1 - m*self.time[t1]
        T = m*t + b
        return T #Temperatures at at that time

    def getTemperatures(self,t,granularity):
        start = self.position[0]
        end = self.position[len(self.position)-1]
        xcoords = np.linspace(start,end,end*granularity)
        T = self.interpolateTempTime(t)
        Tnew = np.interp(xcoords,self.position,T)
        return Tnew

    def getTimeIndices(self,t):
        t1 = 1.0*t/self.tstep
        if(t1 == int(t1)):
            t2 = int(t1)
        else:
            t2 = int(t1 + 1)
        t1 = int(t1)
        return (t1,t2)

    def plotTemperature(self,t,granularity):
        start = self.position[0]
        end = self.position[len(self.position)-1]
        print("start: ", start, " end: ", end, " granularity: ", granularity)
        xcoords = np.linspace(start,end,end*granularity)
        (t1,t2) = self.getTimeIndices(t)
        T1 = self.Temperatures[t1]
        T2 = self.Temperatures[t2]
        T = self.interpolateTempTime(t)
        Tnew = self.getTemperatures(t,granularity)

        plt.plot(self.position,T1,label = ("Time: ",  self.time[t1]))
        plt.plot(self.position,T2,label = ("Time: ",  self.time[t2]))
        plt.plot(self.position,T,label = ("Time: ", t))
        plt.xlabel("Position")
        plt.ylabel("Temperature")
        plt.scatter(xcoords,Tnew)
        plt.legend(loc = 'upper left')
        plt.title("Interpolation Temperature")
        plt.show()
