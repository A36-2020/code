import matplotlib.pyplot as plt
import numpy as np
import math as m
from mpl_toolkits import mplot3d
from F100 import *
from Moment_of_Inertia import *
from material import *
# values
#-------------------------------
Nz = 81
Nx = 41
Ca = 0.505
la = 1.611
Izz_total, Iyy_total = Moment_of_inertia()

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
z = np.zeros((Nz,1))
x = np.zeros((1,Nx))

for i in range (Nz):
    thetaz = (i/Nz)*m.pi
    thetaz2 = ((i+1)/Nz)*m.pi

    zPos = -0.5*((Ca/2)*(1-m.cos(thetaz))+(Ca/2)*(1-m.cos(thetaz2)))
    z[i,0] = zPos


for j in range (Nx):
    thetax = (j/Nx)*m.pi
    thetax2 = ((j+1)/Nx)*m.pi

    xPos = 0.5*((la/2)*(1-m.cos(thetax))+(la/2)*(1-m.cos(thetax2)))
    x[0,j] = xPos

#---------------------------------------


# Integrating aerodynamic load along z-axis and calculating center of pressure
#Returns q integration for every x-axis segment
# simpson integration method
def simpson(F,z,x):
    a = z[0]
    b = z[-1]
    CoP = []
    Q = []
    dz = abs(float((b-a)/2))
    for i in np.arange(Nx):
        q = []
        cop = []
        for j in np.arange(Nz):
            if j % 2 != 0:
                q.append(4*F[j,i])
                cop.append(4*(F[j,i]*z[j]))
            elif j % 2 == 0:
                q.append(2*F[j,i])
                cop.append(2*(F[j,i]*z[j]))
        int_f = (F[0,i]+sum(q)+F[-1,i])*dz/3
        int_fz= (F[0,i]*a+sum(cop)+F[-1,i]*b)*dz/3
        Q.append(int_f)
        CoP.append(int_fz/int_f)
    return Q,CoP

def interpolation_over_span(x_to_interpolate,xcoord_list, q, qz, spanlength):
    #This function interpolates the values for Q and the z-coord of q on a list of xcoords.
    #These coords can range from 0-1, and have to be given in list form.
    #The values for Q are also scaled; meaning that if more points are supplied, the distributed load is discretised
    #in more smaller point loads.

    oringinal_interval_len = spanlength/len(q)
    new_interval_len = spanlength/len(x_to_interpolate)
    qlist = []
    qzlist = []
    for i in x_to_interpolate:
        xcoord  = i*spanlength
        lower_than_x = []
        lower_than_x_qz = []
        higher_than_x = []
        higher_than_x_qz = []
        fractionfound = False
        lower_qcoord = 0
        for j in range(len(q)):
            qcoord = xcoord_list[0][j]
            if qcoord <= xcoord:
                lower_than_x.append(q[j])
                lower_than_x_qz.append(qz[j])
                lower_qcoord = qcoord
            if qcoord > xcoord:
                if fractionfound == False:
                    fraction = (xcoord-lower_qcoord)/(qcoord-lower_qcoord)
                    fractionfound = True
                higher_than_x.append(q[j])
                higher_than_x_qz.append(qz[j])

        lower_qz = qz[0]
        upper_qz = qz[-1]
        lower_q = q[0]
        upper_q = q[-1]

        if len(lower_than_x) > 0:
            lower_qz = lower_than_x_qz[-1]
            lower_q = lower_than_x[-1]
        if len(higher_than_x) > 0:
            upper_qz = higher_than_x_qz[0]
            upper_q = higher_than_x[0]

        interpolated_q = (lower_q+(upper_q-lower_q)*fraction) * new_interval_len/oringinal_interval_len
        interpolated_qz = lower_qz +(upper_qz-lower_qz)*fraction

        qlist.append(interpolated_q)
        qzlist.append(interpolated_qz)
    return qlist, qzlist


interpolated_xlist = np.linspace(0,1,400)
interpolated_xlist = interpolated_xlist[1:]


Qthingy,CoP = simpson(arr,z,x)
print("Qthingy:",np.sum(np.asarray(Qthingy)))
plt.scatter(x[0],Qthingy)
plt.show()
a,b = interpolation_over_span(interpolated_xlist,x,Qthingy,CoP,la)

print("Sum of interpolated data:",np.sum(np.asarray(a)))

scaled_interpolated_xlist = []
for i in interpolated_xlist:
    scaled_interpolated_xlist.append(i*la)

def Mz_Q(x1,x2,x3,Q,x,CoP,Ca):
    Q = (np.asarray(Q))
    x = np.asarray(x)
    CoP = np.asarray(CoP)

    mask_x1 = (x[0]<x1)
    mask_x2 = (x[0]<x2)
    mask_x3 = (x[0]<x3)

    Q_x1 = Q[mask_x1]
    print(Q_x1)
    Q_x2 = Q[mask_x2]
    Q_x3 = Q[mask_x3]

    x_x1 = x[0][mask_x1]
    print(x_x1)
    x_x2 = x[0][mask_x2]
    x_x3 = x[0][mask_x3]

    Mz_Q    = Q*x
    Mz_Q_x1 = Q_x1*(x1-x_x1)**3
    Mz_Q_x2 = Q_x2*(x2-x_x2)**3
    Mz_Q_x3 = Q_x3*(x3-x_x3)**3
    Mx_Q    = Q*(-CoP-0.5*h)

    return Mz_Q, Mz_Q_x1, Mz_Q_x2, Mz_Q_x3, Mx_Q

def moments_x_Q(q,qx,x):
    q = -(np.asarray(q))
    momentsum = 0 #[Nm]
    for i in range(len(qx)):
        qvalue  = -q[i]
        qxvalue = qx[i]
        momentsum = momentsum + qvalue * (qxvalue-x)
    return float(momentsum)


Mz_Q, Mz_Q_x1, Mz_Q_x2, Mz_Q_x3, Mx_Q = Mz_Q(x1,x2,x3,Qthingy,x,Ca,CoP)
#print("Moments:")
#print(Mz_Q, Mz_Q_x1, Mz_Q_x2, Mz_Q_x3, Mx_Q)

def bigmatrix(P,x1,x2,x3,xa,ca,ha,E,Izz_total,Iyy_total,theta,Qthingy,Mx_Q,Mz_Q,Mz_Q_x1,Mz_Q_x2,Mz_Q_x3):

    #Matrix with unknowns (left part)

    bm = np.zeros((11,11))
    #Row 1:
    #Sum of forces in y
    bm[0][0] = 1
    bm[0][1] = 1
    bm[0][2] = 1
    bm[0][6] = m.sin(theta/180.*m.pi)
    #Row 2:
    #Sum of forces in z
    bm[1][3] = 1
    bm[1][4] = 1
    bm[1][5] = 1
    bm[1][6] = m.cos(theta/180.*m.pi)

    #Row 3
    #Sum of moments around x
    bm[2][3] = 0
    bm[2][4] = 0
    bm[2][5] = 0
    bm[2][6] = - m.cos(theta/180.*m.pi)*ha/2 + m.sin(theta/180.*m.pi)*ha/2
    #Row 4
    #Sum of moments around y
    bm[3][3] = x1
    bm[3][4] = x2
    bm[3][5] = x3
    bm[3][6] = m.cos(theta/180.*m.pi)*(x2-xa/2)
    #Row 5
    #Sum of moments around z
    bm[4][0] = x1
    bm[4][1] = x2
    bm[4][2] = x3
    bm[4][6] = m.sin(theta/180.*m.pi) * (x2 - xa / 2)

    #Row 6
    bm[5][7] = x1/(E*Iyy_total)
    bm[5][8] = 1/(E*Iyy_total)

    #Row 7
    bm[6][0] = -(x2-x1)**3/(6*E*Iyy_total)
    bm[6][6] = -m.sin(theta/180.*m.pi)*(xa/2)**3/(6*E*Iyy_total)
    bm[6][7] = 1/(E*Iyy_total)*x2
    bm[6][8] = 1/(E*Iyy_total)
    #Row 8
    bm[7][0] = -(x3-x1)**3/(6*E*Iyy_total)
    bm[7][1] = -(x3-x2)**3/(6*E*Iyy_total)
    bm[7][6] = -m.sin(theta/180.*m.pi)*(x3-(x2-xa/2))**3/(6*E*Iyy_total)
    bm[7][7] = 1/(E*Iyy_total)*x3
    bm[7][8] = 1/(E*Iyy_total)
    #Row 9
    bm[8][9]= x1/(E*Izz_total)
    bm[8][10]= 1/(E*Izz_total)

    #Row 10
    bm[9][3] = (x2-x1)**3/(6*E*Izz_total)
    bm[9][6] = m.cos(theta/180.*m.pi)*(xa/2)**3/(6*E*Izz_total)
    bm[9][9] = x2/(E*Izz_total)
    bm[9][10] = 1/(E*Izz_total)

    #Row 11
    bm[10][3] = (x3-x1)**3/(6*E*Izz_total)
    bm[10][4] = (x3-x2)**3/(6*E*Izz_total)
    bm[10][6]=  m.cos(theta/180.*m.pi)*(x3-(x2-xa/2))**3/(6*E*Izz_total)
    bm[10][9] = x3/(E*Izz_total)
    bm[10][10] = 1/(E*Izz_total)

    #Matrix of knowns (right part)

    bm_knowns = np.zeros((11,1))

    bm_knowns[0] = P*m.sin(theta/180.*m.pi)-np.sum(np.asarray(Qthingy))
    bm_knowns[1] = P*m.cos(theta/180.*m.pi)
    bm_knowns[2] = - P*m.cos(theta/180.*m.pi)*ha/2 + P*m.sin(theta/180.*m.pi)*ha/2 + np.sum(Mx_Q)
    bm_knowns[3] = P*m.cos(theta/180.*m.pi)*(x2+xa/2)
    bm_knowns[4] = P*m.sin(theta/180.*m.pi)*(x2+xa/2)-np.sum(Mz_Q)
    bm_knowns[5] = 1/(6*E*Iyy_total)*np.sum(Mz_Q_x1)+d1
    bm_knowns[6] = np.sum(Mz_Q_x2)/(6*E*Iyy_total)
    bm_knowns[7] = np.sum(Mz_Q_x3)/(6*E*Iyy_total)+P*m.sin(theta/180.*m.pi)*(x3-(x2+xa/2))**3*1/(6*E*Iyy_total)+d3
    bm_knowns[8] = 0
    bm_knowns[9] = 0
    bm_knowns[10]= P*m.cos(theta/180.*m.pi)*(x3-(x2+xa/2))**3/(6*E*Izz_total)

    print("BM:")
    print(bm)
    print("BM_knowns:")
    print(bm_knowns)

    variables = np.linalg.solve(bm,bm_knowns)
    R1y = variables[0]
    R2y = variables[1]
    R3y = variables[2]
    R1z = variables[3]
    R2z = variables[4]
    R3z = variables[5]
    A   = variables[6]
    C1y = variables[7]
    C2y = variables[8]
    C1z = variables[9]
    C2z = variables[10]
    return R1y,R2y,R3y,R1z,R2z,R3z,A,C1y,C2y,C1z,C2z
print("")
print(P,x1,x2,x3,xa,Ca,h,E,Izz_total,Iyy_total,theta,sum(Mx_Q),sum(Mz_Q_x1),sum(Mz_Q_x2),sum(Mz_Q_x3))
print("")
R1y,R2y,R3y,R1z,R2z,R3z,A,C1y,C2y,C1z,C2z = bigmatrix(P,x1,x2,x3,xa,Ca,h,E,Iyy_total,Izz_total,theta,Qthingy,Mx_Q,Mz_Q,Mz_Q_x1,Mz_Q_x2,Mz_Q_x3)
print("SOLUTIONS")
print(R1y,R2y,R3y,R1z,R2z,R3z,A,C1y,C2y,C1z,C2z)

def shear_force_in_y_calculations(R1y, x1, R2y, x2, R3y, x3, A, P, xa, Qvalues, la, xsteps):
    interpolated_xlist=np.linspace(0,la,xsteps)
    scaled_down_interpolated_xlist = np.linspace(0,1,xsteps)
    shear_due_to_aero, b = interpolation_over_span(scaled_down_interpolated_xlist, x, Qvalues, CoP, la)

    shearvalues = [0]
    R1added = False
    R2added = False
    R3added = False
    Padded = False
    Aadded = False
    for i in range(len(interpolated_xlist)):
        localshear = shearvalues[-1] + shear_due_to_aero[i]

        if interpolated_xlist[i]>x1 and R1added == False:
            R1added = True
            localshear = localshear + R1y
        if interpolated_xlist[i]>x2 and R2added == False:
            R2added = True
            localshear = localshear + R2y
        if interpolated_xlist[i]>x3 and R3added == False:
            R3added = True
            localshear = localshear + R3y
        if interpolated_xlist[i]>x2-xa/2 and Aadded == False:
            Aadded = True
            localshear = localshear+A*m.sin(30/180.*m.pi)
        if interpolated_xlist[i]>x2+xa/2 and Padded == False:
            Padded = True
            localshear = localshear-P*m.sin(30/180.*m.pi)

        shearvalues.append(localshear)
    shearvalues = shearvalues[1:]
    return shearvalues,interpolated_xlist

def internal_moment_calculations(shear, shear_x):
    internal_moment_list = [0]
    xdelta = shear_x[1] - shear_x[0]
    for i in range(len(shear_x)):
        internal_moment = internal_moment_list[-1] + shear[i]*xdelta
        internal_moment_list.append(internal_moment)
    internal_moment_list = internal_moment_list[1:]
    return internal_moment_list

def shear_force_in_z_calculations(R1x, x1, R2x, x2, R3x, x3, A, P, xa, la, xsteps):
    interpolated_xlist=np.linspace(0,la,xsteps)
    shearvalues = [0]
    R1added = False
    R2added = False
    R3added = False
    Padded = False
    Aadded = False
    for i in range(len(interpolated_xlist)):
        localshear = shearvalues[-1]
        if interpolated_xlist[i] > x1 and R1added == False:
            R1added = True
            localshear = localshear + R1x
        if interpolated_xlist[i] > x2 and R2added == False:
            R2added = True
            localshear = localshear + R2x
        if interpolated_xlist[i] > x3 and R3added == False:
            R3added = True
            localshear = localshear + R3x
        if interpolated_xlist[i] > x2 - xa / 2 and Aadded == False:
            Aadded = True
            localshear = localshear + A * m.cos(30 / 180. * m.pi)
        if interpolated_xlist[i] > x2 + xa / 2 and Padded == False:
            Padded = True
            localshear = localshear - P * m.cos(30 / 180. * m.pi)

        shearvalues.append(localshear)
    shearvalues = shearvalues[1:]
    return shearvalues, interpolated_xlist

shear_y, xcoords = shear_force_in_y_calculations(R1y,x1,R2y,x2,R3y,x3,A,P,xa,Qthingy,la,200)
shear_z,xcoords = shear_force_in_z_calculations(R1z, x1, R2z, x2, R3z, x3, A, P, xa, la, 200)
internal_moments_in_z = internal_moment_calculations(shear_z, xcoords)
internal_moments_in_y = internal_moment_calculations(shear_y, xcoords)
plt.plot(xcoords,internal_moments_in_z)
plt.show()