#helper functions for R conversion

import numpy as np

def rolling_window(array, window):
    shape = array.shape[:-1] + (array.shape[-1] - window + 1, window)
    strides = array.strides + (array.strides[-1],)
    return np.lib.stride_tricks.as_strided(array, shape=shape, strides=strides)

def rolling_average(array,window):
    return np.mean(rolling_window(array,window),-1)

def rolling_mult(array,window,factor):
    print("len: " , len(array))
    if(len(array) > 1000000):
        full = len(array)
        half = full/2
        window_array1 = rolling_window(array[0,half],window)
        window_array2 = rolling_window(array[half,full])
        tmp1 = window_array1*factor
        tmp2 = window_array2*factor
        tmp = np.concatenate(tmp1,tmp2)
        return tmp
    else:
        return np.sum(rolling_window(array,window)*factor,-1)
