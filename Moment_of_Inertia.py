import math
import numpy as np
from Centroid import *
from F100 import *


def Moment_of_inertia():
    Izz_skin= 2*((tsk * (areaskin/tsk)**3*math.sin(angle)**2)/12 + areaskin*(y_centroid[0]-ycentroidtotal)**2) #for total straight skin
    Iyy_skin= 2*((tsk * (areaskin/tsk)**3*math.cos(angle)**2)/12 + areaskin*(z_centroid[0]-zcentroidtotal)**2) #for total straight skin

    Izz_spar= tsp*h**3/12+ areaspar*(y_centroid[3]-ycentroidtotal)**2
    Iyy_spar= h*tsp**3/12+ areaspar*(z_centroid[3]-zcentroidtotal)**2

    Izz_circ= (math.pi*(0.5*h)**3*tsk)/2 + areacirc*(y_centroid[1]-ycentroidtotal)**2
    #Iyy_circ= (math.pi/8 - 8/(9*math.pi))*((0.5*h)**4-(0.5*h-tsk)**4) + areacirc*(z_centroid[1]-zcentroidtotal)**2
    # Iyy_circ = 0.5*(8+math.pi**2)*tsk*(0.5*h)**3/(4*math.pi) + areacirc*(z_centroid[1]-zcentroidtotal)**2
    Iyy_circ= (math.pi*(0.5*h)**3*tsk)/2 - areacirc*(2*0.5*h/math.pi)**2+ areacirc*(z_centroid[1]-zcentroidtotal)**2

    Izz_total = Izz_skin + Izz_spar + Izz_circ
    Iyy_total = Iyy_skin + Iyy_spar + Iyy_circ

    Izz_stringer = 0
    Iyy_stringer = 0
    for i in range(11):
        Izz_stringer+= areast*(y_centroid[4+i]-ycentroidtotal)**2
        Iyy_stringer+= areast*(z_centroid[4+i]-zcentroidtotal)**2

    Izz_total+= Izz_stringer
    Iyy_total+= Iyy_stringer

    return (Izz_total, Iyy_total)

Izz_total, Iyy_total= Moment_of_inertia()

def Moment_of_Inertia_triangle():
    Izz_skin_t = 2*((tsk * (areaskin/tsk)**3*math.sin(angle)**2)/12 + areaskin*(y_centroids_straight[0])**2) #for total straight skin
    Iyy_skin_t = 2 * ((tsk * (areaskin / tsk) ** 3 * math.cos(angle) ** 2) / 12 + areaskin * (z_centroids_straight[0]-0.5*h)**2)  # for total straight skin

    Izz_spar_t = tsp * h ** 3 / 12 + areaspar * (y_centroids_straight[2]) ** 2
    Iyy_spar_t = h * tsp ** 3 / 12 + areaspar * (z_centroids_straight[2] - 0.5*h) ** 2

    Izz_total_t = Izz_skin_t + Izz_spar_t
    Iyy_total_t = Iyy_skin_t + Iyy_spar_t

    Izz_stringer_t = 0
    Iyy_stringer_t = 0
    for i in range(8):
        Izz_stringer_t += areast * (y_centroids_straight[3 + i]) ** 2
        Iyy_stringer_t += areast * (z_centroids_straight[3 + i] - 0.5*h) ** 2

    Izz_total_t += Izz_stringer_t
    Iyy_total_t += Iyy_stringer_t

    return Izz_total_t, Iyy_total_t


def Moment_of_Inertia_semicirc():
    Izz_circ_c = (math.pi * (0.5 * h) ** 3 * tsk) / 2 + areacirc * (y_centroids_semicirc[0]) ** 2
    Iyy_circ_c = (math.pi / 8 - 8 / (9 * math.pi)) * ((0.5 * h) ** 4 - (0.5 * h - tsk) ** 4) + areacirc * (
                z_centroids_semicirc[0] - 0.5*h) ** 2

    Izz_spar_c = tsp * h ** 3 / 12 + areaspar * (y_centroids_semicirc[1]) ** 2
    Iyy_spar_c = h * tsp ** 3 / 12 + areaspar * (z_centroids_semicirc[1] - 0.5 * h) ** 2

    Izz_total_c = Izz_circ_c + Izz_spar_c
    Iyy_total_c = Iyy_spar_c + Iyy_circ_c

    Izz_stringer_c = 0
    Iyy_stringer_c = 0
    for i in range(3):
        Izz_stringer_c += areast * (y_centroids_straight[2 + i]) ** 2
        Iyy_stringer_c += areast * (z_centroids_straight[2 + i] - 0.5 * h) ** 2

    Izz_total_c += Izz_stringer_c
    Iyy_total_c += Iyy_stringer_c

    return Izz_total_c, Iyy_total_c

def calcJ():
    ## Calculation of J
    dy = h / 2.
    A1 = math.pi * dy ** 2 / 2.
    A2 = (Ca - dy) * dy

    A = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
    b = np.array([0., 0., 0.])

    ## Row 1
    A[0, 0] = 2. * A1
    A[0, 1] = 2. * A2
    b[0] = 1

    ## Row 2
    A[1, 0] = (dy * math.pi / tsk + 2 * dy / tsp) / (2 * A1)
    A[1, 1] = (-2 * dy / tsp) / (2 * A1)
    A[1, 2] = -1.
    b[1] = 0.

    ## Row 3
    A[2, 0] = (-2 * dy / tsp) / (2 * A2)
    A[2, 1] = (2 * areaskin / tsk**2 + 2 * dy / tsp) / (2 * A2)
    A[2, 2] = -1
    b[2] = 0.

    solution = np.linalg.solve(A, b)
    J = 1. / solution[-1]
    return J

if  __name__ == '__main__':
    print(calcJ())
    print(Moment_of_inertia())
    print(Moment_of_Inertia_triangle())
    print(Moment_of_Inertia_semicirc())
