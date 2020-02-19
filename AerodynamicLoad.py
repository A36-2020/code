import matplotlib.pyplot as plt
import numpy as np
import math as m
from mpl_toolkits import mplot3d
from Integrations import *
# values
#-------------------------------
Nz = 81
Nx = 41
Ca = 0.505
la = 1.611

#array with different loads
#------------------------------
lst=[]
lst2=[]
f = open("aerodynamicloadf100.dat","r")
text = f.readlines()
f.close()

for line in text:
    element = line.split(",")
    lst.append(element)

ar = np.array(lst)
arr = ar.astype(np.float)

# connect with correct position
#-------------------------------
z = np.zeros(Nz)
x = np.zeros(Nx)
for i in range (Nz):
    thetaz = (i/Nz)*m.pi
    thetaz2 = ((i+1)/Nz)*m.pi

    zPos = -0.5*((Ca/2)*(1-m.cos(thetaz))+(Ca/2)*(1-m.cos(thetaz2)))
    z[i] = zPos# print(area)



for j in range (Nx):
    thetax = (j/Nx)*m.pi
    thetax2 = ((j+1)/Nx)*m.pi

    xPos = 0.5*((la/2)*(1-m.cos(thetax))+(la/2)*(1-m.cos(thetax2)))
    x[0] = xPos

#Calculate resultant q and CoP for every aileron segment along x-axis (spanwise)
#--------------------------------
area, CoP = simpson(arr,z,x)
