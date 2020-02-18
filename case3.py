import numpy as np
from math import *
import matplotlib.pyplot as plt

import F100
import material
from geometry import A_l, A_t

def case3(x_vals, M_1):
    """ x is a numpy array and M1 is the virtual moment from case 2
        returns internal Moment due to flex on every x as np.array
        as well as the shear flow per section and per cell"""

    # 1 is jammed, number 2 has the force of P
    x_act_1 = F100.x2-F100.xa/2
    x_act_2 = F100.x2+F100.xa/2

    # figure out the moment arm of the actuators taken around the hingeline
    d = F100.h/2 * (1-sin(radians(F100.theta)))
    F_2 = F100.P
    # P is positive Moment, M1 as well, so F_1 = 1/d(M_1+P*d)
    F_1 = M_1/d + F_2

    moment = []
    for x in x_vals:
        v = 0.
        if x >= x_act_1:
            v -= F_1*d
        if x >= F100.x2:
            v += M_1
        if x >= x_act_2:
            v += F_2*d
        moment.append(v)
    moment = np.array(moment)

    # from internal moment we can do angle of twist and shear flow
    S_l = F100.h/F100.tsp+pi*F100.h/F100.tsk
    S_t = F100.h/F100.tsp+2*sqrt(F100.h**2/4+(F100.Ca-F100.h/2)**2)/F100.tsk
    T_l_t = S_t*A_l()**2/S_l/A_t()**2

    q_t = []
    q_l = []

    for m in moment:
        T_t = m/(1+T_l_t)
        T_l = T_l_t*T_t

        q_t.append(T_t/2/A_t())
        q_l.append(T_l/2/A_t())

    q_t = np.array(q_t)
    q_l = np.array(q_l)

    return (moment, q_l, q_t, F_1)
    
l = np.linspace(0, F100.la, 10000)
a = case3(l, 8000)


plt.subplot(221)
plt.plot(l, a[0])
plt.subplot(222)
plt.plot(l, a[1])
plt.subplot(223)
plt.plot(l, a[2])
plt.subplot(224)
plt.plot(l, a[2] - a[1])
plt.show()
plt.show()
plt.show()
plt.show()
