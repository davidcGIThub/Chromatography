#helper functions for R
import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import as_strided as strided

def pyRNorm(n, mean, std):
    nn = int(n)
    arr = np.random.normal(mean,std,len(std))
    np.random.shuffle(arr)
    arr = arr[0:nn]
    return arr

def rolling_window(array, window):
    shape = array.shape[:-1] + (array.shape[-1] - window + 1, window)
    strides = array.strides + (array.strides[-1],)
    return np.lib.stride_tricks.as_strided(array, shape=shape, strides=strides)

def rolling_average(array,window):
    return np.mean(rolling_window(array,window),-1)

def running_average_cumsum(seq, window=100):
    s = np.insert(np.cumsum(seq), 0, [0])
    return (s[window :] - s[:-window]) * (1. / window)

def rolling_mult2(array,factor):
    window = len(factor)
    return np.sum(rolling_window(array,window)*factor,-1)

def rolling_mult(array,factor):
    length = len(array)
    window = len(factor)
    if window > length:
        print ("Error: factor length is larger than array")
    else:
        df = pd.Series(data=array)
        df = df.rolling(window).apply(lambda x: (x*factor).sum() , raw = True)
        df = df.values[window-1:length]
        return df
