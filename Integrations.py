import numpy as np


# simpson interpolation
def simpson(f,z):
    f = -f
    a = z[0]
    b = z[-1]
    n = len(f)
    q = []
    CoP = []
    dx = abs(float((b-a)/2))
    # print(dx)
    for i in np.arange(n):
        if i % 2 != 0:
            q.append(4*f[i])
            CoP.append(4*f[i]*z[i])
        elif i % 2 == 0:
            q.append(2*f[i])
            CoP.append(2*f[i]*z[i])
    int_f = (f[0]+sum(q)+f[-1])*dx/3
    int_fz= (f[0]*a+sum(CoP)+f[-1]*b)*dx/3
    return int_f, (int_fz/int_f)
