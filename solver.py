import matplotlib.pyplot as plt
import numpy as np
import math as m
from mpl_toolkits import mplot3d
from F100 import *
from Moment_of_Inertia import *
from material import *
from integrate import *
from new_shear import section
# values
#-------------------------------
Nz = 81
Nx = 41
Ca = 0.505
la = 1.611
SC = -0.0816
#d1 = 0
#d3 = 0

Izz_total, Iyy_total = Moment_of_inertia()
print(Iyy_total,Izz_total)
J = calcJ()


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

    original_interval_len = spanlength/len(q)
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

        interpolated_q = (lower_q+(upper_q-lower_q)*fraction) * new_interval_len/original_interval_len
        interpolated_qz = lower_qz +(upper_qz-lower_qz)*fraction

        qlist.append(-interpolated_q)
        qzlist.append(interpolated_qz)
    return qlist, qzlist

interpolated_xlist_len = 1000
interpolated_xlist  = np.zeros(interpolated_xlist_len)
for j in range (interpolated_xlist_len):
    thetax = (j/interpolated_xlist_len)*m.pi
    thetax2 = ((j+1)/interpolated_xlist_len)*m.pi

    xPos = 0.5*((la/2)*(1-m.cos(thetax))+(la/2)*(1-m.cos(thetax2)))
    interpolated_xlist[j] = xPos/la




Qthingy,CoP = simpson(arr,z,x)
print(sum(Qthingy))
interpolated_qvalues,interpolated_CoPs = interpolation_over_span(interpolated_xlist,x,Qthingy,CoP,la)
print(sum(interpolated_qvalues))
interpolated_xlist = [interpolated_xlist]

#for i in range(len(interpolated_qvalues)):
#    interpolated_qvalues[i] = 0

scaled_interpolated_xlist = []
for i in interpolated_xlist:
    scaled_interpolated_xlist.append(i*la)

def q_y(x,p):
    dx = la/len(Qthingy)
    s = 0
    r = 0
    i = 0
    while s <= (x-dx):
        r += Qthingy[i]*(x-s)**p
        s += dx
        i += 1
    return r#/m.factorial(p)

print(q_y(la,0))

q_y = np.vectorize(q_y)
print(q_y(la,1),q_y(la,0)*(la/2))

def q_t(x,p):
    dx = la/len(Qthingy)
    s = 0
    r = 0
    i = 0
    while s <= (x-dx):
        r += Qthingy[i]*(CoP[i]-SC)
        s += dx
        i += 1
    return r/m.factorial(p)

q_t = np.vectorize(q_t)

def Mz_Q(x1,x2,x3,xA,Q,x,CoP,Ca):
    Q = (np.asarray(Q))
    x = np.asarray(x)
    CoP = np.asarray(CoP)

    mask_x1 = (x[0]<x1)
    mask_x2 = (x[0]<x2)
    mask_x3 = (x[0]<x3)
    mask_xA = (x[0]<xA)

    Q_x1 = Q[mask_x1]
    Q_x2 = Q[mask_x2]
    Q_x3 = Q[mask_x3]
    Q_xA = Q[mask_xA]

    x_x1 = x[0][mask_x1]
    x_x2 = x[0][mask_x2]
    x_x3 = x[0][mask_x3]
    x_xa = x[0][mask_xA]

    Mz_Q    = Q*x
    Mz_Q_x1 = Q_x1*(x1-x_x1)**3
    Mz_Q_x2 = Q_x2*(x2-x_x2)**3
    Mz_Q_x3 = Q_x3*(x3-x_x3)**3
    T_A = Q_xA*(-abs(CoP-SC))
    Mx_A = Q_xA*(xA-x_xa)**3
    Mx_Q    = Q*(-CoP-0.5*h)
    return Mz_Q, Mz_Q_x1, Mz_Q_x2, Mz_Q_x3, Mx_Q,T_A,Mx_A

def moments_x_Q(q,qx,x):
    q = -(np.asarray(q))
    momentsum = 0 #[Nm]
    for i in range(len(qx)):
        qvalue  = -q[i]
        qxvalue = qx[i]
        momentsum = momentsum + qvalue * (qxvalue-x)
    return float(momentsum)

xA = x2-xa/2
Mz_Q, Mz_Q_x1, Mz_Q_x2, Mz_Q_x3, Mx_Q,T_A,Mx_A = Mz_Q(x1,x2,x3,xA,interpolated_qvalues,scaled_interpolated_xlist,Ca,interpolated_CoPs)

#print("Moments:")
#print(Mz_Q, Mz_Q_x1, Mz_Q_x2, Mz_Q_x3, Mx_Q)

def bigmatrix(P,x1,x2,x3,xa,ca,ha,E,Izz_total,Iyy_total,theta,Qthingy,Mx_Q,Mz_Q,Mz_Q_x1,Mz_Q_x2,Mz_Q_x3,G,J,T_A):

    #Matrix with unknowns (left part)
    bm = np.zeros((12,12))

    bm_knowns = np.zeros((12, 1))

    #Row 1:
    #Sum of forces in y
    bm[0][0] = 1
    bm[0][1] = 1
    bm[0][2] = 1
    bm[0][6] = m.sin(theta/180.*m.pi)

    bm_knowns[0] = P * m.sin(theta / 180. * m.pi) - np.sum(np.asarray(Qthingy))

    #Row 2:
    #Sum of forces in z
    bm[1][3] = 1
    bm[1][4] = 1
    bm[1][5] = 1
    bm[1][6] = m.cos(theta/180.*m.pi)

    bm_knowns[1] = P * m.cos(theta / 180. * m.pi)

    #Row 3
    #Sum of moments around x
    bm[2][3] = 0
    bm[2][4] = 0
    bm[2][5] = 0
    bm[2][6] = - m.cos(theta/180.*m.pi)*ha/2 + m.sin(theta/180.*m.pi)*ha/2

    bm_knowns[2] = -P * m.cos(theta / 180. * m.pi) * ha / 2 + P * m.sin(theta / 180. * m.pi) * ha / 2 + np.sum(Mx_Q)

    #Row 4
    #Sum of moments around y
    bm[3][3] = x1
    bm[3][4] = x2
    bm[3][5] = x3
    bm[3][6] = m.cos(theta/180.*m.pi)*(x2-xa/2)

    bm_knowns[3] = P * m.cos(theta / 180. * m.pi) * (x2 + xa / 2)

    #Row 5
    #Sum of moments around z
    bm[4][0] = x1
    bm[4][1] = x2
    bm[4][2] = x3
    bm[4][6] = m.sin(theta/180.*m.pi) * (x2 - xa / 2)

    bm_knowns[4] = P * m.sin(theta / 180. * m.pi) * (x2 + xa / 2) + np.sum(Mz_Q)

    #Row 6
    bm[5][7] = x1
    bm[5][8] = 1

    bm_knowns[5] = 1 / (6 * E * Iyy_total) * np.sum(Mz_Q_x1) + d1 * m.cos(theta / 180. * m.pi)

    #Row 7
    bm[6][0] = -(x2-x1)**3/(6*E*Iyy_total)
    bm[6][6] = -m.sin(theta/180.*m.pi)*(xa/2)**3/(6*E*Iyy_total)
    bm[6][7] = x2
    bm[6][8] = 1

    bm_knowns[6] = np.sum(Mz_Q_x2) / (6 * E * Iyy_total)

    #Row 8
    bm[7][0] = -(x3-x1)**3/(6*E*Iyy_total)
    bm[7][1] = -(x3-x2)**3/(6*E*Iyy_total)
    bm[7][6] = -m.sin(theta/180.*m.pi)*(x3-(x2-xa/2))**3/(6*E*Iyy_total)
    bm[7][7] = x3
    bm[7][8] = 1

    bm_knowns[7] = np.sum(Mz_Q_x3) / (6 * E * Iyy_total) + P * m.sin(theta / 180. * m.pi) * (
                x3 - (x2 + xa / 2)) ** 3 * 1 / (6 * E * Iyy_total) + d3 * m.cos(theta / 180. * m.pi)

    #Row 9
    bm[8][9]= x1
    bm[8][10]= 1

    bm_knowns[8] = -d1 * m.sin(theta / 180. * m.pi)

    #Row 10
    bm[9][3] = (x2-x1)**3/(6*E*Izz_total)
    bm[9][6] = m.cos(theta/180.*m.pi)*(xa/2)**3/(6*E*Izz_total)
    bm[9][9] = x2
    bm[9][10] = 1

    bm_knowns[9] = 0
#
    #Row 11
    bm[10][3] = (x3-x1)**3/(6*E*Izz_total)
    bm[10][4] = (x3-x2)**3/(6*E*Izz_total)
    bm[10][6]=  m.cos(theta/180.*m.pi)*(x3-(x2-xa/2))**3/(6*E*Izz_total)
    bm[10][9] = x3
    bm[10][10] = 1

    bm_knowns[10] = P * m.cos(theta / 180. * m.pi) * (x3 - (x2 + xa / 2)) ** 3 / (6 * E * Izz_total) - d3 * m.sin(
        theta / 180. * m.pi)

    #Row 12
    bm[11][0] = 1/6*(E*Iyy_total)*((x2-xa/2)-x1)**3*m.sin(theta/180.*m.pi)
    bm[11][4] = SC*m.sin(theta/180.*m.pi)+1/(6*E*Izz_total)*((x2-xa/2)-x1)**3*m.cos(theta/180.*m.pi)
    bm[11][11] = 1

    bm_knowns[11] = (np.sum(T_A) * SC - 1 / (6 * E * Iyy_total) * np.sum(Mx_A)) * m.sin(theta / 180. * m.pi)

    #print(bm)

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
    CT = variables[11]
    return R1y,R2y,R3y,R1z,R2z,R3z,A,C1y,C2y,C1z,C2z,CT


vals = bigmatrix(P,x1,x2,x3,xa,Ca,h,E,Iyy_total,Izz_total,theta,interpolated_qvalues,Mx_Q,Mz_Q,Mz_Q_x1,Mz_Q_x2,Mz_Q_x3,G,J,T_A)
print(vals)

def maucaly(x, xn):
    return np.where(x>xn,x-xn,0)

def maucaly0(x, xn):
    return np.where(x>xn,1,0)

def deflectiony(x, one, two, three, A, P, C1, C2):
    return -1/6/E/Izz_total*(one*maucaly(x,x1)**3+A*m.sin(theta/180*m.pi)*maucaly(x,x2-xa/2)**3+two*maucaly(x,x2)**3-P*m.sin(theta/180*m.pi)*maucaly(x,x2+xa/2)**3 +three*maucaly(x,x3)**3-q_y(x,3))+C1*x+C2

def sheary(x, one, two, three, A, P, C1, C2):
    return one*maucaly0(x,x1)+A*m.sin(theta/180*m.pi)*maucaly0(x,x2-xa/2)+two*maucaly0(x,x2)-P*m.sin(theta/180*m.pi)*maucaly0(x,x2+xa/2) +three*maucaly0(x,x3)-q_y(x,0)

def momenty(x, one, two, three, A, P, C1, C2):
    return one*maucaly(x,x1)**1+A*m.sin(theta/180*m.pi)*maucaly(x,x2-xa/2)**1+two*maucaly(x,x2)**1-P*m.sin(theta/180*m.pi)*maucaly(x,x2+xa/2)**1 +three*maucaly(x,x3)**1-q_y(x,1)


def deflectionz(x, one, two, three, A, P, C1, C2):
    return -1/6/E/Iyy_total*(one*maucaly(x,x1)**3+A*m.cos(theta/180*m.pi)*maucaly(x,x2-xa/2)**3+two*maucaly(x,x2)**3-P*m.cos(theta/180*m.pi)*maucaly(x,x2+xa/2)**3 +three*maucaly(x,x3)**3)+C1*x+C2

def shearz(x, one, two, three, A, P, C1, C2):
    return +one*maucaly0(x,x1)+A*m.cos(theta/180*m.pi)*maucaly0(x,x2-xa/2)+two*maucaly0(x,x2)-P*m.cos(theta/180*m.pi)*maucaly0(x,x2+xa/2) +three*maucaly0(x,x3)

def momentz(x, one, two, three, A, P, C1, C2):
    return +one*maucaly(x,x1)**1-A*m.cos(theta/180*m.pi)*maucaly(x,x2-xa/2)**1+two*maucaly(x,x2)**1+P*m.cos(theta/180*m.pi)*maucaly(x,x2+xa/2)**1 +three*maucaly(x,x3)**1
a  = np.linspace(0,la,1000)

def twist(x, oney, twoy, threy, A, P, C):
    s = -SC
    t = theta/180*m.pi
    return 1/G/J*(-oney*s*maucaly(x,x1)-twoy*s*maucaly(x,x2)-threy*s*maucaly(x,x3)-P*(m.sin(t)*s-m.cos(t)*h/2)*maucaly(x,x2+xa/2)+A*(m.cos(t)*h/2-m.sin(t)*s)*maucaly(x,x2-xa/2)-twoy*s*maucaly(x,x3)+q_t(x,0))

#-47000, 65000, -18000
plt.plot(a, sheary(a,vals[0],vals[1],vals[2], vals[6], P, vals[7], vals[8]))
plt.ylabel("Shear Force y [N]")
plt.xlabel("X position [m]")
plt.show()
plt.plot(a, momenty(a, vals[0],vals[1],vals[2], vals[6], P, vals[7], vals[8]))
plt.ylabel("Internal Moment z [Nm]")
plt.xlabel("X position [m]")
plt.show()
plt.plot(a, deflectiony(a, vals[0],vals[1],vals[2], vals[6], P, vals[7], vals[8]))
plt.ylabel("Deflection y [m]")
plt.xlabel("X position [m]")
plt.show()

plt.plot(a, shearz(a,vals[3],vals[4],vals[5], vals[6], P, vals[9], vals[10]))
plt.ylabel("Shear Force z [N]")
plt.xlabel("X position [m]")
plt.show()
plt.plot(a, momentz(a,vals[3],vals[4],vals[5], vals[6], P, vals[9], vals[10]))
plt.ylabel("Internal Moment y [Nm]")
plt.xlabel("X position [m]")
plt.show()
plt.plot(a, deflectionz(a,vals[3],vals[4],vals[5], vals[6], P, vals[9], vals[10]))
plt.ylabel("Deflection z [m]")
plt.xlabel("X position [m]")
plt.show()

plt.plot(a, twist(a,vals[0],vals[1],vals[2], vals[6], P, vals[11]))
plt.show()

x = 0.45
Vy =  sheary(x,vals[0],vals[1],vals[2], vals[6], P, vals[7], vals[8])[0]
My =  momenty(x,vals[0],vals[1],vals[2], vals[6], P, vals[7], vals[8])[0]
Vz = shearz(x,vals[3],vals[4],vals[5], vals[6], P, vals[9], vals[10])[0]
print("dasfasdf ", Vz)
Mz = momentz(x,vals[3],vals[4],vals[5], vals[6], P, vals[9], vals[10])[0]
T = 1185
s = section(0.00005, Vy,Vz,T,My,Mz)
s.show()
