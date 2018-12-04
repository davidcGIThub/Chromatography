import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

res = abs(np.loadtxt("resolution.csv",delimiter = ','))
pos = np.loadtxt("pos_x.csv",delimiter = ',')
vel = np.loadtxt("vel_x.csv",delimiter = ',')
std = np.loadtxt("std.csv",delimiter = ',')
molecExit = np.loadtxt("molecExitTime.csv",delimiter = ',')
t = np.loadtxt("timeOutput.csv",delimiter = ',')

plt.figure(1)
for i in range(len(pos)):
    plt.plot(t,pos[i],label = "Compound " + str(i))
plt.legend(loc = 'best')
plt.xlabel("Time (sec)")
plt.ylabel("Molecule Location")
plt.title("Average Position of Compounds")
plt.show()

plt.figure(2)
for i in range(len(vel)):
    plt.plot(t,vel[i],label = "Compound " + str(i))
plt.legend(loc = 'best')
plt.xlabel("Time (sec)")
plt.ylabel("Molecule Velocity")
plt.title("Average Velocity of Compounds")
plt.show()

plt.figure(3)
for i in range(len(std)):
    plt.plot(t,std[i],label = "Compound " + str(i))
plt.legend(loc = 'best')
plt.xlabel("Time (sec)")
plt.ylabel("Molecule Standard Deviation")
plt.title("Standard Deviation of Each Compound")
plt.show()

plt.figure(4)
if len(np.shape(res))==1:
    plt.plot(t,res,label = "Compound " + str(i) + " to " + str(i+1))
else:
    for i in range(len(res)):
        plt.plot(t,res[i],label = "Compound " + str(i) + " to " + str(i+1))
plt.legend(loc = 'best')
plt.xlabel("Time (sec)")
plt.ylabel("Resolution")
plt.title("Resolution")
plt.show()

plt.figure(5)
for i in range(len(molecExit)):
    y = np.linspace(0,1,len(molecExit[0]))
    print("len y: ", len(y))
    print("molecExit[i]: ", len(molecExit[i]))
    np.random.shuffle(y)
    plt.scatter(molecExit[i],y,label = "Compound " + str(i))
plt.legend(loc = 'best')
plt.xlabel("Time (sec)")
plt.ylabel("Exit Time")
plt.title("Compound Exit Time")
plt.show()

plt.figure(6)
for i in range(len(molecExit)):
    sns.kdeplot(molecExit[i],color='k')
plt.xlim(0,max(molecExit.flatten())+5)
plt.yticks([])
plt.show()
