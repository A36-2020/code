import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Wedge
from F100 import *

z_centroid= []
y_centroid= []
Area=[]

angle = math.atan((0.5*h)/(Ca-0.5*h))
stspace = (2*math.sqrt((0.5*h)**2+ (Ca-0.5*h)**2)+ math.pi*0.5*h)/11
zcircfront = stspace*(0.5*h)/(0.5*math.pi*0.5*h)
y_centroidcirc = math.sqrt((0.5*h)**2-((0.5*h)- zcircfront)**2) #this is for the stringers on the circular part
z_centroidcirc = zcircfront

z_centroidspar = 0.5*h
y_centroidspar = 0

z_centroidskinup = Ca- (Ca-0.5*h)/2
y_centroidskinup = 0.25*h

z_centroidskindown = Ca- (Ca-0.5*h)/2
y_centroidskindown = -0.25*h

z_centroid_semicirc = 0.5*h - 2*0.5*h/math.pi
y_centroid_semicirc = 0

y_centroid+=[y_centroidskinup, y_centroid_semicirc, y_centroidskindown, y_centroidspar]
z_centroid+=[z_centroidskinup, z_centroid_semicirc, z_centroidskindown, z_centroidspar]
      
areaskin = math.sqrt((0.5*h)**2+(Ca-0.5*h)**2)*tsk
areacirc = math.pi*0.5*h*tsk
areaspar = h*tsp
areast = tst*hst + tst*wst

Area+=[areaskin, areacirc, areaskin, areaspar]
Area+= 11*[areast]


def stringercentroids(y_centroid, z_centroid):
    for n in range(4):
        y_centroidpos = math.sin(angle)*(0.5*stspace+n*stspace)
        y_centroid.append(y_centroidpos)
        z_centroidstraight = Ca- math.cos(angle)*(0.5*stspace+(n)*stspace)
        z_centroid.append(z_centroidstraight)
        
    z_centroid+=[z_centroidcirc, 0, z_centroidcirc]
    y_centroid+=[y_centroidcirc, 0, -y_centroidcirc]

    for n in range(4):
        y_centroidneg = math.sin(angle)*(0.5*stspace+(3-n)*stspace)
        y_centroid.append(-y_centroidneg)
        z_centroidstraight = Ca- math.cos(angle)*(0.5*stspace+(3-n)*stspace)
        z_centroid.append(z_centroidstraight)
        
    return y_centroid, z_centroid

#print(stringercentroids(y_centroid, z_centroid))
#print(Area)

ztimesA = 0
ytimesA = 0

y_centroid, z_centroid = stringercentroids(y_centroid, z_centroid)
for i in range(len(Area)):
    ztimesA += z_centroid[i]*Area[i]
    ytimesA += y_centroid[i]*Area[i]

areasum = sum(Area)

zcentroidtotal = ztimesA/areasum
ycentroidtotal = ytimesA/areasum

print(z_centroid)
if __name__ == '__main__':
    print(zcentroidtotal, ycentroidtotal)

zup = [0.5*h,Ca]
yup = [0.5*h,0]
zdown = [0.5*h,Ca]
ydown = [-0.5*h,0]
zspar = [0.5*h,0.5*h]
yspar = [0.5*h, -0.5*h]



plt.scatter(zcentroidtotal, ycentroidtotal, label="Centroid", color= "blue", marker= "x", s=35)
plt.scatter(z_centroid, y_centroid, label="Stiffeners", color = "orange", marker= "x", s=25)
plt.scatter(z_centroid_semicirc, y_centroid_semicirc, label= "Semicircle", color= "green", marker= "x", s=30)
plt.scatter(z_centroidspar, y_centroidspar, label= "Spar", color= "red", marker= "x", s=30)
plt.scatter(z_centroidskinup, y_centroidskinup, label= "Skins", color= "purple", marker= "x", s=30)
plt.scatter(z_centroidskindown, y_centroidskindown, color= "purple",  
            marker= "x", s=30)
plt.plot(zup,yup, color="black", linewidth=0.7)
plt.plot(zdown,ydown, color="black", linewidth=0.7)
plt.plot(zspar, yspar, color="black", linewidth=0.1)
circle1 = Wedge((0.5*h,0), 0.5*h, 90, 270, color="black", fill = False)
fig = plt.gcf()
ax = fig.gca()
ax.add_artist(circle1)
fig.set_size_inches(8,3)

plt.legend()
plt.show()
