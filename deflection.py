import numpy as np
import matplotlib.pyplot as plt

del1 = 0.389/100
del3 = 1.245/100
l1 = 0.498-0.125
l2 = 1.494-0.498

# ax^3 + bx^2 + cx
# dx^3 + bx^2 + cx

A = np.zeros((4,4))

A[0,:] = [(-l1)**3, + (-l1)**2, -l1, 0]
A[1,:] = [0, l2**2, l2, l2**3]
A[2,:] = [6*(-l1), 2*(-l1), 0, 0]
A[3,:] = [0, 2*l2, 0, 6*l2]

b = [del1, del3, 0, 0]

x = np.linalg.solve(A,b)

spline1 = lambda p: x[0]*p**3+x[1]*p**2+x[2]*p 
spline2 = lambda p: x[3]*p**3+x[1]*p**2+x[2]*p 

uno = spline1(np.linspace(0,-l1,50))
dos = spline2(np.linspace(l2,0,50))

plt.plot(np.linspace(0,-l1,50), uno)
plt.plot(np.linspace(l2,0,50), dos)
plt.show()
