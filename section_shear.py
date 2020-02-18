import numpy as np
from math import *
import F100
import geometry
import material
import matplotlib.pyplot as plt

class section():
    def __init__(self, n, V):
        self.n = n
        I = geometry.Inertia_xx()
        V_lt = geometry.S_al()/geometry.S_at()
        V_t = V/(1+V_lt)
        V_l = V_lt*V_t
        print(V_l, V_t)

        self.p_l = np.linspace(0,geometry.S_l(),n)
        self.p_t = np.linspace(0,geometry.S_t(),n)

        self.q_l = np.array([self.line_integral_l(x) for x in self.p_l])
        self.q_t = np.array([self.line_integral_t(x) for x in self.p_t])

        q_l0 = self.integrate(self.q_l, geometry.S_l()/n)/geometry.S_l()
        q_t0 = self.integrate(self.q_t, geometry.S_t()/n)/geometry.S_t()

        self.q_l += q_l0+q_t0
        self.q_t += q_t0+q_l0
        print(self.q_t[-1], self.q_l[0])

    def show(self):
        x = np.array([self.x_t(p) for p in self.p_t ])
        y = np.array([self.y_t(p) for p in self.p_t ])
        plt.scatter(x,y,c=(self.q_t))
        plt.colorbar()
        plt.axis('equal')
        plt.show()

    def x_l(self, s):
        if s <= F100.h:
            return 0
        
        s -= F100.h
        r = F100.h/2
        phi = s/r
        return -1*sin(phi)*r
    def x_t(self, s):
        if s <= geometry.l():
            return (s*(F100.Ca-F100.h/2)/geometry.l())

        s -= geometry.l()
        if s <= geometry.l():
            return (F100.Ca-F100.h/2)-(s*(F100.Ca-F100.h/2)/geometry.l())
        s -= geometry.l()

        return 0


    def integrate(self, ys, dx):
        return sum([y*dx for y in ys])

    def y_l(self, s):
        if s <= F100.h:
            return s-F100.h/2
        
        s -= F100.h
        r = F100.h/2
        phi = s/r
        return cos(phi)*r


    def line_integral_l(self, end, ds=0.0001):
        t = lambda s: F100.tsp if s <= F100.h else F100.tsk


        integral = 0.0
        s = 0.0
        while s < end:
            integral += ds*t(s)*self.y_l(s)
            s += ds

        return integral

    def y_t(self, s):
        if s <= geometry.l():
            return -F100.h/2+s*(F100.h/2/geometry.l())
        
        s -= geometry.l()
        if s <= geometry.l():
            return s*(F100.h/2/geometry.l())
        s -= geometry.l()

        return -1*(s-F100.h/2)


    def line_integral_t(self, end, ds=0.0001):
        t = lambda s: F100.tsk if s <= 2*geometry.l() else F100.tsp
            
        integral = 0.0
        s = 0.0
        while s < end:
            integral += ds*t(s)*self.y_t(s)
            s += ds

        return integral

s = section(1000, 3000)
s.show()
