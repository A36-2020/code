import numpy as np

def integrate(y,x):
    if y.shape != x.shape:
        return 0
    val = 0.0
    for i in range(1,x.shape[0]):
        val += (y[i-1]+y[i])/2*(x[i]-x[i-1])
    return val

def integrate_2(y,dx):
    val = 0.0
    for i in range(1,y.shape[0]):
        val += (y[i-1]+y[i])/2*(dx)
    return val
