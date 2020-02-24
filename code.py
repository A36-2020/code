import numpy as np


# simpson interpolation
def simpson(f,a,b,n):
    y = []
    dx = (b-a)/2
    for i in np.arange(1,n):
        if i % 2 != 0:
            y.append(4*f[n])
        elif i % 2 == 0:
            y.append(2*f[n])

    return (f[0]+sum(y)+f[-1])*dx/3
