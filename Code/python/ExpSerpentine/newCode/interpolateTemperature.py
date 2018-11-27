import numpy as np
import matplotlib.pyplot as plt

class GradientData:
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

    def getMolecTemp(self,molPos,t):
        molTemp = np.array([])
        for i in range(len(molPos)):
            molTemp = np.append(molTemp,self.getTemperature(molPos[i],t))
        return molTemp

    def getGradient(self,t):
        TempVec = np.array([])
        for i in range(len(self.x)):
            TempVec = np.append(TempVec,self.getTemperature(self.x[i],t))
        return TempVec

    def plotTemperature(self,x,t):
        (x1,x2) = self.getPositionIndices(x)
        (t1,t2) = self.getTimeIndices(t)
        T1 = self.interpolate(self.x[x1],self.x[x2],self.T[t1][x1],self.T[t1][x2],x)
        T2 = self.interpolate(self.x[x1],self.x[x2],self.T[t2][x1],self.T[t2][x2],x)
        T = self.interpolate(self.t[t1],self.t[t2],T1,T2,t)
        plt.plot(self.x,self.T[t1],label = ("Time: ",  self.t[t1]))
        plt.plot(self.x,self.T[t2],label = ("Time: ", self.t[t2]))
        #plt.plot(self.x,self.getGradient(t),label = ("Time: ", t))
        plt.xlabel("Position")
        plt.ylabel("Temperature")
        plt.scatter(x,T,label = ("x:",x," T:",round(T,2),"t:",t))
        plt.legend(loc = 'upper left')
        plt.title("Interpolation Temperature")
        plt.show()
        return T
