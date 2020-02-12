import numpy as np
import matplotlib.pyplot as plt

import F100
def deformation(x_vals):
    """ x is a numpy array, returns deflection due to flex on every x as np.array
        as well as the curvature on every position"""
    l1 = F100.x2 - F100.x1
    l2 = F100.x3 - F100.x2

    # ax^3 + bx^2 + cx
    # dx^3 + bx^2 + cx

    A = np.zeros((4,4))

    A[0,:] = [(-l1)**3, + (-l1)**2, -l1, 0]
    A[1,:] = [0, l2**2, l2, l2**3]
    A[2,:] = [6*(-l1), 2*(-l1), 0, 0]
    A[3,:] = [0, 2*l2, 0, 6*l2]

    b = [F100.d1, F100.d3, 0, 0]

    x = np.linalg.solve(A,b)

    spline1 = lambda p: x[0]*(p-F100.x2)**3+x[1]*(p-F100.x2)**2+x[2]*(p-F100.x2)
    spline2 = lambda p: x[3]*(p-F100.x2)**3+x[1]*(p-F100.x2)**2+x[2]*(p-F100.x2) 

    theta1 = -1*(spline1(F100.x1)-spline1(F100.x1+0.0001))/0.0001
    theta3 = -1*(spline2(F100.x3)-spline2(F100.x3+0.0001))/0.0001

    def deformation(u):
        if u < F100.x1:
            return spline1(F100.x1)+theta1*(u-F100.x1)
        if F100.x1 <= u <= F100.x2:
            return spline1(u)
        if F100.x2 <= u <= F100.x3:
            return spline2(u)
        else:
            return spline2(F100.x3)+theta3*(u-F100.x3)

    defs = []
    for i in x_vals:
        defs.append(deformation(i))

    max_curv = 2*x[1]
    # we know curvature is max at Hinge 2, so and it looks like a tipi is 0 after 

    return np.array(defs)

plt.plot(np.linspace(0, F100.la, 100), deformation(np.linspace(0, F100.la, 100)))
plt.show()
