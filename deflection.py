import numpy as np
import matplotlib.pyplot as plt

import F100


def deformation(x_vals):
    """ x is a numpy array, returns deflection due to flex on every x as np.array
        as well as the Moment times stifness and shear force on every position"""
    l1 = F100.x2 - F100.x1
    l2 = F100.x3 - F100.x2

    # ax^3 + bx^2 + cx
    # dx^3 + bx^2 + cx

    A = np.zeros((4, 4))

    A[0, :] = [(-l1) ** 3, ((-l1) ** 2), -l1, 0]
    A[1, :] = [0, l2 ** 2, l2, l2 ** 3]
    A[2, :] = [6 * (-l1), 2, 0, 0]
    A[3, :] = [0, 2, 0, 6 * l2]

    b = [F100.d1, F100.d3, 0, 0]

    x = np.linalg.solve(A, b)

    spline1 = (
        lambda p: x[0] * (p - F100.x2) ** 3
        + x[1] * (p - F100.x2) ** 2
        + x[2] * (p - F100.x2)
    )
    spline2 = (
        lambda p: x[3] * (p - F100.x2) ** 3
        + x[1] * (p - F100.x2) ** 2
        + x[2] * (p - F100.x2)
    )

    theta1 = -1 * (spline1(F100.x1) - spline1(F100.x1 + 0.0001)) / 0.0001
    theta3 = -1 * (spline2(F100.x3) - spline2(F100.x3 + 0.0001)) / 0.0001

    def deformation(u):
        if u < F100.x1:
            return spline1(F100.x1) + theta1 * (u - F100.x1)
        if F100.x1 <= u <= F100.x2:
            return spline1(u)
        if F100.x2 <= u <= F100.x3:
            return spline2(u)
        else:
            return spline2(F100.x3) + theta3 * (u - F100.x3)

    defs = []
    for i in x_vals:
        defs.append(deformation(i))

    defs = np.array(defs)

    return (defs, np.gradient(np.gradient(defs, x_vals), x_vals) ,np.gradient(np.gradient(np.gradient(defs, x_vals), x_vals), x_vals))

l = np.linspace(0, F100.la, 10000)
a = deformation(l)
plt.subplot(212)
plt.plot(l, a[0])

plt.subplot(221)
plt.plot(l, a[1])

plt.subplot(222)
plt.plot(l, a[2])

plt.show()
