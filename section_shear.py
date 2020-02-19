import numpy as np
from math import *
import F100
import Moment_of_Inertia
import geometry
import material
import matplotlib.pyplot as plt

class section():
    def __init__(self, ds, V, V_z, T):
        self.d = ds
        I = Moment_of_Inertia.Moment_of_inertia()[0]
        V_lt = geometry.S_al()/geometry.S_at()
        V_t = V/(1+V_lt)
        V_l = V_lt*V_t
        V_t_z = V_z/(1+V_lt)
        V_l_z = V_lt*V_t_z
        print(V_l, V_t, V_l + V_t, V)

        self.p_l = np.arange(0,geometry.S_l(),ds)
        self.p_t = np.arange(0,geometry.S_t(),ds)

        self.q_l = V_l*np.array([self.line_integral_l(x) for x in self.p_l])/Moment_of_Inertia.Moment_of_Inertia_semicirc()[0]+V_l_z*np.array([self.line_integral_l_z(x) for x in self.p_l])/Moment_of_Inertia.Moment_of_Inertia_semicirc()[1]
        self.q_t = V_t*np.array([self.line_integral_t(x) for x in self.p_t])/Moment_of_Inertia.Moment_of_Inertia_triangle()[1]+V_t_z*np.array([self.line_integral_t_z(x) for x in self.p_t])/Moment_of_Inertia.Moment_of_Inertia_triangle()[1]

        q_l0 = self.qs0_l()
        q_t0 = self.qs0_t()

        self.q_l += q_l0
        self.q_t += q_t0

        #DO the tourqe
        T_l_t = geometry.S_at()*geometry.A_l()**2/geometry.S_al()/geometry.A_t()**2

        T_t = T/(1+T_l_t)
        T_l = T_l_t*T_t

        self.q_t += T_t/2/geometry.A_t()
        self.q_l += T_l/2/geometry.A_l()



        # extract the spar
        l_end_i = int(F100.h/ds)
        self.s = -1*(self.q_l[:l_end_i]-(self.q_t[-l_end_i:]))

        self.p_x = np.array([0 for _ in self.s])
        self.p_y = np.linspace(-F100.h/2,F100.h/2,self.s.shape[0])

        self.q_l = self.q_l[l_end_i:]
        self.p_l = self.p_l[l_end_i:]

        self.q_t = self.q_t[:-l_end_i]
        self.p_t = self.p_t[:-l_end_i]

        print(self.s[0],self.s[-1],self.q_l[0],self.q_l[-1],self.q_t[0],self.q_t[-1])
        input()
        
        
    def show(self):
        xt = np.array([self.x_t(p) for p in self.p_t ])
        yt = np.array([self.y_t(p) for p in self.p_t ])
        xl = np.array([self.x_l(p) for p in self.p_l ])
        yl = np.array([self.y_l(p) for p in self.p_l ])

        x = np.hstack((xt,xl,self.p_x))
        y = np.hstack((yt,yl,self.p_y))
        c = np.hstack((self.q_t, self.q_l, self.s))

        m = max(np.abs(c))
        plt.scatter(x,y,c=(c))#, vmin=-m, vmax=m)
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

    def line_integral_l_z(self, end, ds=0.0001):
        t = lambda s: F100.tsp if s <= F100.h else F100.tsk


        integral = 0.0
        s = 0.0
        while s < end:
            integral += ds*t(s)*self.x_l(s)
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

    def line_integral_t_z(self, end, ds=0.0001):
        t = lambda s: F100.tsk if s <= 2*geometry.l() else F100.tsp
            
        integral = 0.0
        s = 0.0
        while s < end:
            integral += ds*t(s)*self.x_t(s)
            s += ds

        return integral

    def qs0_l(self):
        t = lambda s: F100.tsp if s <= F100.h else F100.tsk

        integral = 0.0
        integral2 = 0.0
        for i in range(self.q_l.shape[0]):
            integral += self.d * self.q_l[i]/t(self.p_l[i])
            integral2 += self.d/t(self.p_l[i])
        return integral/integral2


    def qs0_t(self):
        t = lambda s: F100.tsk if s <= 2*geometry.l() else F100.tsp

        integral = 0.0
        integral2 = 0.0
        for i in range(self.q_t.shape[0]):
            integral += self.d * self.q_t[i]/t(self.p_t[i])
            integral2 += self.d/t(self.p_t[i])
        return integral/integral2

s = section(0.001, 000,3000, 00)
s.show()
