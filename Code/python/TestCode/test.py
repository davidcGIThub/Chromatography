#test file
import HelperFunctionsR as hfr
import matplotlib.pyplot as plt
import numpy as np
import random as rnd
import seaborn as sns
import pandas as pd
from ctypes import cdll
import time
j = 2.0 #number of compounds
n_mol = 100.0
#
# t = time.time()
# gradient = pd.read_csv("ColGrad0001s.csv",header=None).values
# Temp_vec = gradient[:,1].flatten()
# xpos = gradient[:,0].flatten()
# print('Loading from pandas took {}s'.format(time.time() - t))
#
# t = time.time()
# index, npvalues = np.load('ColGrad0001s.npz')
# print('Loading from npy took {}s'.format(time.time() - t))
start = time.time()
f = pd.DataFrame()
for i in range(1000):
  f['col_{}'.format(i)] = np.random.rand(25000)
print('Generating data took {}s'.format(time.time() - start))

start = time.time()
f.to_hdf('frame.hdf', 'main', format='fixed')
print('Saving to HDF took {}s'.format(time.time() - start))

start = time.time()
np.savez('frame.npz', f.index, f.values)
print('Saving to npy took {}s'.format(time.time() - start))

start = time.time()
pd.read_hdf('frame.hdf')
print('Loading from HDF took {}s'.format(time.time() - start))

start = time.time()
index, values = np.load('frame.npz')
print('Loading from npy took {}s'.format(time.time() - start))
print(values)
