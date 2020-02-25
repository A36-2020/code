from math import *

from AerodynamicLoad import *
from F100 import *

<<<<<<< HEAD

=======
print(len(Q),len(x[0]))
>>>>>>> wieger
mask_x1 = (x[0]<x1)
mask_x2 = (x[0]<x2)
mask_x3 = (x[0]<x3)

Q_x1 = Q[mask_x1]
Q_x2 = Q[mask_x2]
Q_x3 = Q[mask_x3]

x_x1 = x[0][mask_x1]
x_x2 = x[0][mask_x2]
x_x3 = x[0][mask_x3]

Mz_Q_x1 = Q_x1*x_x1
Mz_Q_x2 = Q_x2*x_x2
Mz_Q_x3 = Q_x3*x_x3
