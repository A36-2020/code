import math
import numpy as np
from Centroid import *
from F100 import *

Izz_skin= 2*((tsk * (areaskin/tsk)**3*math.sin(angle)**2)/12 + areaskin*(y_centroid[0]-ycentroidtotal)**2) #for total straight skin
Iyy_skin= 2*((tsk * (areaskin/tsk)**3*math.cos(angle)**2)/12 + areaskin*(z_centroid[0]-zcentroidtotal)**2) #for total straight skin

Izz_spar= tsp*h**3/12+ areaspar*(y_centroid[3]-ycentroidtotal)**2
Iyy_spar= h*tsp**3/12+ areaspar*(z_centroid[3]-zcentroidtotal)**2

Izz_total = Izz_skin + Izz_spar
Iyy_total = Iyy_skin + Iyy_spar

Izz_stringer = 0
Iyy_stringer = 0
for i in range(11):
    Izz_stringer+= areast*(y_centroid[4+i]-ycentroidtotal)**2
    Iyy_stringer+= areast*(z_centroid[4+i]-zcentroidtotal)**2

Izz_total+= Izz_stringer
Iyy_total+= Iyy_stringer

print(Izz_total, Iyy_total)









