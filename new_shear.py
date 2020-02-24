import numpy as np
import F100
import Moment_of_Inertia
from math import *
import matplotlib.pyplot as plt
from integrate import integrate

c = F100.Ca - F100.h/2
r = F100.h/2
l = sqrt(r**2+c**2)


def _t(s):
    if s <= 2*l+r*pi:
        return F100.tsk
    else:
        return F100.tsp

t = np.vectorize(_t)

def _arm(s):
    if s <= 2*l:
        return r*c/l
    elif s <= 2*l+r*pi:
        return r
    else:
        return 0

arm = np.vectorize(_arm)

def _line_yy(s):
    val = 0
    if s <= l:
        return F100.tsk*(-r*s+r/2/l*s**2)
    else:
        val += F100.tsk*(-r*l+r/2*l)

    s -= l

    if s <= l:
        return val + F100.tsk*(r/2/l*s**2)
    else:
        val += F100.tsk*(r/2*l)

    s -= l

    if s <= (pi*r):
        return val + F100.tsk * r*r*sin(s/r)
    else:
        val += 0

    s -= pi*r

    return val + F100.tsp*(-r*s+1/2*s**2)

line_yy = np.vectorize(_line_yy)

def _line_zz(s):
    val = 0
    if s <= l:
        return F100.tsk*(c/2/l*s**2)
    else:
        val += F100.tsk*(c/2/l*l**2)

    s -= l

    if s <= l:
        return val + F100.tsk*(c*s-c/2/l*s**2)
    else:
        val += F100.tsk*(c*l-c/2*l)

    s -= l

    if s <= (pi*r):
        return val + F100.tsk * r*r*(cos(s/r)-1)
    else:
        return val + F100.tsk * r*r*(cos(pi)-1)

line_zz = np.vectorize(_line_zz)


def _z(s):
    val = 0
    if s <= l:
        return s*c/l
    else:
        val = c

    s -= l

    if s <= l:
        return val - s*c/l
    else:
        val -= c

    s -= l

    if s <= (pi*r):
        return val - sin(s/r)*r
    else:
        return val

z = np.vectorize(_z)

def _y(s):
    val = 0
    if s <= l:
        return -r+s*r/l
    else:
        val = 0

    s -= l

    if s <= l:
        return val + s*r/l
    else:
        val += r

    s -= l

    if s <= (pi*r):
        return cos(s/r)*r
    else:
        val = -r

    s -= pi*r

    return val + s

y = np.vectorize(_y)


Iyy, Izz = Moment_of_Inertia.Moment_of_inertia()

class section():
    def __init__(self, ds, Vy, Vz, T, M=0):

        endt = int(2*l/ds)
        endl = endt + int(r*pi/ds)

        self.s = np.arange(0,2*l+r*pi+2*r,ds)
        self.qsy = Vy*line_yy(self.s)/Iyy
        self.qsz = Vz*line_zz(self.s)/Izz
        self.qs = self.qsz + self.qsy

        self.qb = -1*integrate((self.qs/t(self.s)), self.s)/integrate((1/t(self.s)), self.s)

        self.z = z(self.s)
        self.y = y(self.s)




        q = self.qs+self.qb
        
        print(integrate(q*arm(self.s)/t(self.s),self.s))

        m = np.max(np.abs(q))
        
        plt.scatter(self.z,self.y,c=q, cmap='seismic', vmin=-m, vmax=m)
        plt.axis('equal')
        plt.colorbar()
        plt.show()

section(0.001, 1, 0, 0)
