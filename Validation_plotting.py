from nbconvert.exporters import notebook
import xlrd
from mpl_toolkits import mplot3d
import numpy as np
import math as m
import matplotlib.tri as mtri
import scipy as sp
import pylab as pl
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from matplotlib import cm

""" Importing the Excel file where all the validation data has been formatted and saved"""

loc = ("Validation_data.xlsx")

wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)

"""Saving the x, y and z coordinates of all the nodes"""

X_nodes = np.array([float(sheet.cell_value(0, 1))])
Y_nodes = np.array([float(sheet.cell_value(0, 2))])
Z_nodes = np.array([float(sheet.cell_value(0, 3))])

for i in range(sheet.nrows-1):
    X_nodes = np.append(X_nodes, float(sheet.cell_value(i + 1, 1)))
    Y_nodes = np.append(Y_nodes, float(sheet.cell_value(i + 1, 2)))
    Z_nodes = np.append(Z_nodes, float(sheet.cell_value(i + 1, 3)))

"""Creating an array with all the x y z coordinates for each node, without deflection"""

all_nodes = np.array([[X_nodes[0], Y_nodes[0], Z_nodes[0]]])

for i in range(len(X_nodes)-1):
    all_nodes= np.append(all_nodes, [[X_nodes[i+1], Y_nodes[i+1], Z_nodes[i+1]]], axis=0)


"""Creating an array with all x y z coordinates for each node, including deflection"""

sheet_deflection= wb.sheet_by_index(7)

all_nodes_deflections = np.array([[X_nodes[0]+float(sheet_deflection.cell_value(1,2)), Y_nodes[0]+float(sheet_deflection.cell_value(1,3)), Z_nodes[0]+float(sheet_deflection.cell_value(1,4))]])

for i in range(len(X_nodes)-1):
    all_nodes_deflections= np.append(all_nodes_deflections, [[X_nodes[i+1]+float(sheet_deflection.cell_value(i+2,2)), Y_nodes[i+1]+float(sheet_deflection.cell_value(i+2,3)), Z_nodes[i+1]+float(sheet_deflection.cell_value(i+2,4))]], axis=0)


"""Creating an array that stores the numbers of the nodes that generate each element"""

sheet1= wb.sheet_by_index(1)

all_elements = np.array([[int(sheet1.cell_value(1,1)), int(sheet1.cell_value(1,2)), int(sheet1.cell_value(1,3)), int(sheet1.cell_value(1,4))]])

for i in range(sheet1.nrows-2):
    all_elements = np.append(all_elements, [[int(sheet1.cell_value(i+2,1)), int(sheet1.cell_value(i+2,2)), int(sheet1.cell_value(i+2,3)), int(sheet1.cell_value(i+2,4))]], axis=0)


"""Creating an array that stores the coordinates of each of the 4 nodes that generate each element, without deflection"""

all_elements_nodes_coord = np.array([[all_nodes[int(sheet1.cell_value(1,1))-1], all_nodes[int(sheet1.cell_value(1,2))-1], all_nodes[int(sheet1.cell_value(1,3))-1], all_nodes[int(sheet1.cell_value(1,4))-1]]])

for i in range(sheet1.nrows-2):
    all_elements_nodes_coord = np.append(all_elements_nodes_coord, [[all_nodes[int(sheet1.cell_value(i+2, 1)) - 1], all_nodes[int(sheet1.cell_value(i+2, 2)) - 1],
                             all_nodes[int(sheet1.cell_value(i+1, 3)) - 1],
                                all_nodes[int(sheet1.cell_value(i+1, 4)) - 1]]], axis=0)



"""Creating an array that stores the coordinates of each of the 4 nodes that generate each element, with deflection"""

all_elements_nodes_coord_deflect = np.array([[all_nodes_deflections[int(sheet1.cell_value(1,1))-1], all_nodes_deflections[int(sheet1.cell_value(1,2))-1], all_nodes_deflections[int(sheet1.cell_value(1,3))-1], all_nodes_deflections[int(sheet1.cell_value(1,4))-1]]])

for i in range(sheet1.nrows-2):
    all_elements_nodes_coord_deflect = np.append(all_elements_nodes_coord_deflect, [[all_nodes_deflections[int(sheet1.cell_value(i+2, 1)) - 1], all_nodes_deflections[int(sheet1.cell_value(i+2, 2)) - 1],
                             all_nodes_deflections[int(sheet1.cell_value(i+1, 3)) - 1],
                                all_nodes_deflections[int(sheet1.cell_value(i+1, 4)) - 1]]], axis=0)



"""Saving the number of the elements with the highest stress"""
max_stress_spar= 3842
max_stress_skin= 2391
max_stress_overall= 2391

max_shear_spar= 3842
max_shear_skin= 2390
max_shear_overall= 2390


"""Saving all the Von Mises stresses for the elements"""

sheet_stress= wb.sheet_by_index(6)
von_misses_stresses= np.array([float(sheet_stress.cell_value(1,6))])
for i in range(sheet_stress.nrows-2):
    von_misses_stresses= np.append(von_misses_stresses, [float(sheet_stress.cell_value(i+2,6))], axis=0 )



"""Saving all the shear stresses for the elements"""

shear_stresses= np.array([float(sheet_stress.cell_value(1,7))])
for i in range(sheet_stress.nrows-2):
    shear_stresses= np.append(shear_stresses, [float(sheet_stress.cell_value(i+2,7))], axis=0 )


"""Plotting all the nodes in 3D"""

def plotting_all_nodes():

    fig = plt.figure()
    ax = plt.axes(projection="3d")

    ax.scatter(X_nodes, Y_nodes, Z_nodes, c='r', marker='o')
    ax.set_title('All nodes')


    plt.show()


#vtx = sp.rand(4,3)
#print(vtx)
#print(all_elements_nodes_coord_deflect[0])

"""Defining a function for normalisation"""

def normal(x, x_min, x_max):
    return (x-x_min)/(x_max-x_min)

"""Plotting the von misses stress distribution along the elements in 3D, by plotting each element as a separate polygon"""

def plot_polygon_von_misses():
    ax = a3.Axes3D(pl.figure())

    for i in range(len(all_elements_nodes_coord_deflect)):

        vtx1 = all_elements_nodes_coord_deflect[i,0]
        vtx2 = all_elements_nodes_coord_deflect[i,1]
        vtx3 = all_elements_nodes_coord_deflect[i,2]
        vtx4 = all_elements_nodes_coord_deflect[i,3]
        if max(abs(vtx1-vtx2))<100 and  max(abs(vtx2-vtx3))<100 and max(abs(vtx3-vtx4))<100 and max(abs(vtx2-vtx4))<100 and max(abs(vtx3-vtx1))<100 and max(abs(vtx4-vtx1))<100:
            vtx = [vtx1,vtx2,vtx4,vtx3]
            tri = a3.art3d.Poly3DCollection([vtx])
            col = normal(von_misses_stresses[i], min(von_misses_stresses), max(von_misses_stresses))
            cmap=plt.get_cmap("plasma")
            col=[cmap(col)]
            tri.set_color(col)
            ax.add_collection3d(tri)


    ax.set(xlim=(0,3000), ylim=(-150,150), zlim=(-600, 200))
    pl.show()

"""Plotting the shear stress distribution along the elements in 3D, by plotting each element as a separate polygon"""
def plot_polygon_shear():
    fig = pl.figure()
    ax = a3.Axes3D(fig)
    m = cm.ScalarMappable(cmap=plt.get_cmap("plasma"))
    for i in range(len(all_elements_nodes_coord_deflect)):

        vtx1 = all_elements_nodes_coord_deflect[i,0]
        vtx2 = all_elements_nodes_coord_deflect[i,1]
        vtx3 = all_elements_nodes_coord_deflect[i,2]
        vtx4 = all_elements_nodes_coord_deflect[i,3]
        if max(abs(vtx1-vtx2))<100 and  max(abs(vtx2-vtx3))<100 and max(abs(vtx3-vtx4))<100 and max(abs(vtx2-vtx4))<100 and max(abs(vtx3-vtx1))<100 and max(abs(vtx4-vtx1))<100:
            vtx = [vtx1,vtx2,vtx3,vtx4]
            tri = a3.art3d.Poly3DCollection([vtx])
            col = normal(shear_stresses[i], min(shear_stresses), max(shear_stresses))
            cmap=plt.get_cmap("plasma")
            col=[cmap(col)]
            tri.set_color(col)
            ax.add_collection3d(tri)
    m.set_array(np.array([min(shear_stresses), max(shear_stresses)]))
    ax.set(xlim=(0,3000), ylim=(-150,150), zlim=(-600, 200))
    fig.colorbar(m)

    pl.show()

#plot_polygon_von_misses()
#plot_polygon_shear()

"""Saving an "average" point for each element -> the coordinates of the point that would be exactly in the middle of each rectangular element"""

point_for_element= [(all_elements_nodes_coord[0][0]+all_elements_nodes_coord[0][1]+all_elements_nodes_coord[0][2]+all_elements_nodes_coord[0][3])/4]
for i in range(len(all_elements_nodes_coord)-1):
    point_for_element.append((all_elements_nodes_coord[i+1][0]+all_elements_nodes_coord[i+1][1]+all_elements_nodes_coord[i+1][2]+all_elements_nodes_coord[i+1][3])/4)

point_for_element=np.array(point_for_element)


"""Plotting the von mises stress distribution for the cross-section at the position of the element with highest von mises stress"""

def cross_section_max_von_misses(x_coord):
    fig = plt.gcf()
    ax = fig.gca()
    ax.invert_xaxis()
    cmap = plt.get_cmap("jet")
    m = cm.ScalarMappable(cmap=plt.get_cmap("jet"))
    col=[]
    nodes_list=[]
    c=0
    for node in list(point_for_element):
        if abs(node[0]-x_coord)<15 and (abs(node[1]+node[2])>5 or abs(node[1]-node[2])==0) and int(node[1])!=46 and abs(node[1]+46.277)>2:
            nodes_list.append(list(node))
            col1 = normal(von_misses_stresses[c], min(von_misses_stresses), max(von_misses_stresses))

            col.append(cmap(col1))
        c+=1

    nodes_array= np.array(nodes_list)



    for i in range(len(nodes_array)):

        plt.scatter(nodes_array[i][2], nodes_array[i][1], c=col[i], marker='o')

    m.set_array(np.array([min(von_misses_stresses), max(von_misses_stresses)]))
    fig.colorbar(m)
    ax.set_title('Von Mises stress distribution')
    ax.set_ylabel('y (m)')
    ax.set_xlabel('z (m)')
    plt.show()

"""Plotting the shear stress distribution for the cross-section at the position of the element with highest shear stress"""

def cross_section_max_shear_stress(x_coord):
    fig = plt.gcf()
    ax = fig.gca()
    ax.invert_xaxis()
    cmap = plt.get_cmap("jet")
    m = cm.ScalarMappable(cmap=plt.get_cmap("jet"))
    col=[]
    nodes_list=[]
    c=0
    for node in list(point_for_element):
        if abs(node[0]-x_coord)<15 and (abs(node[1]+node[2])>5 or abs(node[1]-node[2])==0) and int(node[1])!=46 and abs(node[1]+46.277)>2:
            nodes_list.append(list(node))
            col1 = normal(shear_stresses[c], min(shear_stresses), max(shear_stresses))

            col.append(cmap(col1))
        c+=1

    nodes_array= np.array(nodes_list)



    for i in range(len(nodes_array)):

        plt.scatter(nodes_array[i][2], nodes_array[i][1], c=col[i], marker='o')

    m.set_array(np.array([min(shear_stresses), max(shear_stresses)]))
    fig.colorbar(m)
    ax.set_title('Shear stress distribution')
    ax.set_ylabel('y (m)')
    ax.set_xlabel('z (m)')
    plt.show()

cross_section_max_von_misses(point_for_element[2390][0]) #the element with max von mises stress is 2391

#cross_section_max_shear_stress(point_for_element[2389][0]) #the element with max shear stress is 2389

"""Plotting the deflection of the hingeline with respect to y"""

def plotting_deflection_of_hingeline_y():
    nodes_on_hingeline_y=[]
    nodes_on_hingeline_x=[]
    for i in range(len(all_nodes)):
        if all_nodes[i][1]==0 and all_nodes[i][2]==0:
            nodes_on_hingeline_y.append(all_nodes_deflections[i][1])
            nodes_on_hingeline_x.append(all_nodes_deflections[i][0])


    nodes_on_hingeline_y=np.array(nodes_on_hingeline_y)
    nodes_on_hingeline_x = np.array(nodes_on_hingeline_x)

    fig = plt.gcf()
    ax = fig.gca()

    for i in range(len(nodes_on_hingeline_x)):

        plt.scatter(nodes_on_hingeline_x[i], nodes_on_hingeline_y[i], c="b", marker='o')

    ax.set_ylabel('y (mm)')
    ax.set_xlabel('x (mm)')
    ax.set_title('Deflection of hinge line in y')

    plt.show()

"""Plotting the deflection of the hingeline with respect to z"""

def plotting_deflection_of_hingeline_z():
    nodes_on_hingeline_z = []
    nodes_on_hingeline_x = []
    for i in range(len(all_nodes)):
        if all_nodes[i][1] == 0 and all_nodes[i][2] == 0:
            nodes_on_hingeline_z.append(all_nodes_deflections[i][2])
            nodes_on_hingeline_x.append(all_nodes_deflections[i][0])

    nodes_on_hingeline_z = np.array(nodes_on_hingeline_z)
    nodes_on_hingeline_x = np.array(nodes_on_hingeline_x)

    fig = plt.gcf()
    ax = fig.gca()

    for i in range(len(nodes_on_hingeline_x)):
        plt.scatter(nodes_on_hingeline_x[i], nodes_on_hingeline_z[i], c="b", marker='o')

    ax.set_ylabel('z (mm)')
    ax.set_xlabel('x (mm)')
    ax.set_title('Deflection of hinge line in z')

    plt.show()

mini=1000
c=100000
for i in range(len(all_nodes)):

    if int(all_nodes_deflections[i][1])==int(all_nodes_deflections[i][2]) and int(all_nodes_deflections[i][1])==0:
        if all_nodes_deflections[i][1]<mini or all_nodes_deflections[i][2]<mini:
            mini=min(all_nodes_deflections[i][1], all_nodes_deflections[i][2])
            c=i

#print(all_nodes_deflections[c], c)
min_deflec=all_nodes[c]

"""Saving nodes on the leading edge"""
nodes_on_LE_z=[]
nodes_on_LE_y=[]
nodes_LE=[]
for i in range(len(all_nodes)):
    if all_nodes[i][1] == 0 and all_nodes[i][2] == max(Z_nodes):
        nodes_LE.append(i)
        nodes_on_LE_y.append(all_nodes_deflections[i][1])
        nodes_on_LE_z.append(all_nodes_deflections[i][2])

"""Saving nodes on the trailing edge"""
nodes_on_TE_z=[]
nodes_on_TE_y=[]
nodes_TE=[]
for i in range(len(all_nodes)):
    if all_nodes[i][1] == 0 and all_nodes[i][2] == min(Z_nodes):
        nodes_TE.append(i)
        nodes_on_TE_y.append(all_nodes_deflections[i][1])
        nodes_on_TE_z.append(all_nodes_deflections[i][2])

"""Saving pairs of points LE-TE to draw chord and extract twist angle"""
pairs_LE_TE=[]

for i in range(len(nodes_LE)):
    for j in range(len(nodes_TE)):
        if all_nodes[nodes_LE[i]][0]==all_nodes[nodes_TE[j]][0]:
            pairs_LE_TE.append([nodes_LE[i], nodes_TE[j]])


"""Plotting the twist angle"""
def plot_twist():
    twist_along_x=[]
    x_coords=[]
    for i in range(len(pairs_LE_TE)):
        twist_along_x.append(m.atan(all_nodes_deflections[pairs_LE_TE[i][0]][1]/all_nodes_deflections[pairs_LE_TE[i][0]][2]))
        x_coords.append(all_nodes[pairs_LE_TE[i][0]][0])

    fig = plt.gcf()
    ax = fig.gca()

    for i in range(len(pairs_LE_TE)):
        plt.scatter(x_coords[i], twist_along_x[i], c="b", marker='o')

    ax.set_ylabel('twist')
    ax.set_xlabel('x (mm)')
    ax.set_title('Twist along the hinge line')

    plt.show()

#plot_twist()
#plotting_deflection_of_hingeline_y()

#def plotting_twist():

