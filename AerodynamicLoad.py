import matplotlib.pyplot as plt
import numpy as np
import math as m
from mpl_toolkits import mplot3d

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
f = open("Aileronload.dat","r")
text = f.readlines()
f.close()

for line in text:
    element = line.split(",")
    lst.append(element)

ar = np.array(lst)
arr = ar.astype(np.float)


# connect with correct position
#-------------------------------
z = np.zeros((80,1))
x = np.zeros((1,40))

for i in range (80):
    thetaz = (i/Nz)*m.pi
    thetaz2 = ((i+1)/Nz)*m.pi

    zPos = -0.5*((Ca/2)*(1-m.cos(thetaz))+(Ca/2)*(1-m.cos(thetaz2)))
    z[i,0] = zPos


for j in range (40):
    thetax = (j/Nx)*m.pi
    thetax2 = ((j+1)/Nx)*m.pi

    xPos = 0.5*((la/2)*(1-m.cos(thetax))+(la/2)*(1-m.cos(thetax2)))
    x[0,j] = xPos

# plot surface
#---------------------------------------



#output array sample:
q = [1,2,5,4,2,3]
qz = [0.2,0.2,0.4,0.6,0.2,0.5]
xlist = [0.1,0.2,0.4,0.6,0.8,1]
x = [0.25,0.5,0.75,1]
x = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
spanlength = 12

def interpolation_over_span(x, q, qz, spanlength):
    oringinal_interval_len = spanlength/len(q)
    new_interval_len = spanlength/len(x)
    qlist = []
    qzlist = []
    for i in x:
        xcoord  = i*spanlength
        lower_than_x = []
        lower_than_x_qz = []
        higher_than_x = []
        higher_than_x_qz = []
        fractionfound = False
        for j in range(len(q)):
            qcoord = j/len(q)*spanlength
            if qcoord < xcoord:
                lower_than_x.append(q[j])
                lower_than_x_qz.append(qz[j])
            if qcoord > xcoord:
                if fractionfound == False:
                    fraction = qcoord/xcoord - 1
                    fractionfound = True
                higher_than_x.append(q[j])
                higher_than_x_qz.append(qz[j])

        lower_qz = qz[0]
        upper_qz = qz[-1]
        lower_q = q[0]
        upper_q = q[-1]

        #print(xcoord)

        if len(lower_than_x) > 0:
            lower_qz = lower_than_x_qz[-1]
            lower_q = lower_than_x[-1]
        if len(higher_than_x) > 0:
            upper_qz = higher_than_x_qz[0]
            upper_q = higher_than_x[0]

        #print(lower_qz)
        #print(upper_qz)
        #print(lower_q)
        #print(upper_q)
        #print(fraction)
        #print("")

        interpolated_q = (lower_q+(upper_q-lower_q)*fraction) * new_interval_len/oringinal_interval_len
        interpolated_qz = lower_qz +(upper_qz-lower_qz)*fraction

        #print("Interpolated values q & qz:")
        #print(interpolated_q)
        #print(interpolated_qz)
        #print("")
        #print("")

        qlist.append(interpolated_q)
        qzlist.append(interpolated_qz)
    return qlist, qzlist

a,b = interpolation_over_span(x,q,qz,spanlength)
print(a)
print(b)

plt.plot(q,xlist)
plt.scatter(a,x)
plt.show()

r1x = 0.2
r2x = 0.4
r3x = 0.6
E = 1000
I = 1.2

#def maccaulay_reactionforces_inx(r1x,r2x,r3x,x,qlist,spanlength,E,I):






