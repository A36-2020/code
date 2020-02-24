import numpy as np
import F100
import Moment_of_Inertia
from math import *
import matplotlib.pyplot as plt
from integrate import integrate, integrate_2

c = F100.Ca - F100.h/2
r = F100.h/2
l = sqrt(r**2+c**2)


def _t(s):
    if s <= 2*l+r*pi:
        return F100.tsk
    else:
        return F100.tsp

t = np.vectorize(_t)

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
        return val+F100.tsk*(-r*s+r/2/l*s**2)
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
        return val + F100.tsk*(c/2/l*s**2)
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

def _sumy(s):
    val = 0
    if s < l:
        n = floor(s/(l/5))
        y = (1+(n-1)/2)*((l/5))*r/l
        y -= r
        return val + F100.Ast * y * n
    else:
        val = 4*F100.Ast*(-r/2)

    s -= l
    if s < l:
        n = floor(s/(l/5))
        y = (1+(n-1)/2)*((l/5))*r/l
        return val + F100.Ast * y * n
    else:
        val -= 4*F100.Ast*(-r/2)
    s -= l

    if s < r*pi:
        theta = s/r
        if theta > pi/4:
            val += F100.Ast*sin(pi/4)*r
        if theta > 3*pi/4:
            val -= F100.Ast*sin(pi/4)*r

    return val
        
sumy = np.vectorize(_sumy)

def _sumz(s):
    val = 0
    if s < l:
        n = floor(s/(l/5))
        y = (1+(n-1)/2)*((l/5))*c/l
        return val + F100.Ast * y * n
    else:
        val += 4*F100.Ast*(c/2)

    s -= l
    if s < l:
        n = floor(s/(l/5))
        y = (1+(n-1)/2)*((l/5))*c/l
        y = c-y
        return val + F100.Ast * y * n
    else:
        val += 4*F100.Ast*(c/2)
    s -= l

    theta = s/r
    if theta > pi/4:
        val -= F100.Ast*cos(pi/4)*r
    if theta > pi/2:
        val -= F100.Ast*r
    if theta > 3*pi/4:
        val -= F100.Ast*cos(pi/4)*r

    return val
        
sumz = np.vectorize(_sumz)

Iyy, Izz = Moment_of_Inertia.Moment_of_inertia()
zbar = 0.204
print(Iyy, Izz)

class section():
    def __init__(self, ds, Vy, Vz, T, My=0, Mz=0):

        endt = int(2*l/ds)+1
        endl = endt + int(r*pi/ds)

        self.s = np.arange(0,2*l+r*pi+2*r,ds)
        self.qsy = Vy*(line_yy(self.s)+sumy(self.s))/Iyy
        self.qsz = Vz*(line_zz(self.s)+sumz(self.s))/Izz
        self.qs = self.qsz + self.qsy

        qt = self.qs[:endt]
        ql = self.qs[endt:endl]
        sp = self.qs[endl:]

        qbt = -1*(integrate_2(qt, ds)/F100.tsk+integrate_2(-1*sp, ds)/F100.tsp)/(2*l/F100.tsk+2*r/F100.tsp)
        qbl = -1*(integrate_2(ql, ds)/F100.tsk+integrate_2(sp, ds)/F100.tsp)/(r*pi/F100.tsk+2*r/F100.tsp)

        qt += qbt
        ql += qbl
        sp = sp + qbl - qbt

        # Shear flow over, now torque
        A = np.array([[1,1],[0,0]])
        A[1,0] = 4/pi/pi/(r**4)*(r*pi/F100.tsk+2*r/F100.tsp)
        A[1,1] = -1/r/r/c/c*(2*l/F100.tsk+2*r/F100.tsp)

        T1, T2 = np.linalg.solve(A, [T, 0])

        q1 = T1/2/(0.5*r**2*pi)
        q2 = T2/2/(r*c)

        qt += q2
        ql += q1
        sp = sp + q1 - q2


        q = np.hstack((qt,ql,sp))
        # Direct Stress
        sigy = My*y(self.s)/Iyy
        sigz = Mz*(z(self.s)-zbar)/Izz
        self.sig = sigy+sigz

        # Visualisation
        self.z = z(self.s)
        self.y = y(self.s)


        m = np.max(np.abs(q))
        
        # Sanity check bro
        print((q[endt-1]-q[endt]+q[-1])/m)
        print((-q[0]+q[endl-1]-q[endl])/m)

        self.q = q
        self.m = m

        self.mises = np.sqrt(np.square(self.q)+3*np.square(self.q/t(self.s)))

    def show(self):
        plt.scatter(self.z, self.y, c=self.q, vmin=-self.m, vmax=self.m)
        plt.colorbar()
        plt.axis('equal')
        plt.show()
        plt.scatter(self.z, self.y, c=self.sig)
        plt.colorbar()
        plt.axis('equal')
        plt.show()
        plt.scatter(self.z, self.y, c=self.mises)
        plt.colorbar()
        plt.axis('equal')
        plt.show()

s = section(0.0001, 10, 4, 5, 0, 1)
s.show()
