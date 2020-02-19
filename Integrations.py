import numpy as np


# simpson interpolation
def simpson(F,z,x):
    F = -1*F
    a = z[0]
    b = z[-1]
    nz = len(z)
    nx = len(x)
    CoP = []
    Q = []
    dz = abs(float((b-a)/2))
    print(F.shape)
    print(F[0,:])
    for i in np.arange(nx):
        q = []
        cop = []
        for j in np.arange(nz):
            if j % 2 != 0:
                q.append(4*F[j,i])
                cop.append(4*(F[j,i]*z[j]))
            elif j % 2 == 0:
                q.append(2*F[j,i])
                cop.append(2*(F[j,i]*z[j]))
        int_f = (F[0,i]+sum(q)+F[-1,i])*dz/3
        int_fz= (F[0,i]*a+sum(cop)+F[-1,i]*b)*dz/3
        Q.append(int_f)
        CoP.append(int_fz/int_f)
    return Q,CoP
