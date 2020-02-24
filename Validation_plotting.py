from nbconvert.exporters import notebook
import xlrd
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.tri as mtri
import scipy as sp
import pylab as pl
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from matplotlib import cm

loc = ("Validation_data.xlsx")

wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)

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



"""Saving the number of the element with the highest stress"""
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

print(von_misses_stresses[2390])

"""Saving all the shear stresses for the elements"""
shear_stresses= np.array([float(sheet_stress.cell_value(1,7))])
for i in range(sheet_stress.nrows-2):
    shear_stresses= np.append(shear_stresses, [float(sheet_stress.cell_value(i+2,7))], axis=0 )


def plotting_all_nodes():

    fig = plt.figure()
    ax = plt.axes(projection="3d")

    ax.scatter(X_nodes, Y_nodes, Z_nodes, c='r', marker='o')
    ax.set_title('All nodes')


    plt.show()


#vtx = sp.rand(4,3)
#print(vtx)
#print(all_elements_nodes_coord_deflect[0])

def normal(x, x_min, x_max):
    return (x-x_min)/(x_max-x_min)


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
    m.set_array([min(shear_stresses), max(shear_stresses)])
    ax.set(xlim=(0,3000), ylim=(-150,150), zlim=(-600, 200))
    fig.colorbar(m)

    pl.show()

#plot_polygon_von_misses()
plot_polygon_shear()


def nodes_cross_section(x_coord):
    nodes_list=[]
    for node in list(all_nodes):
        if node[0]==x_coord:
            nodes_list.append(list(node))

    nodes_array= np.array(nodes_list)
    print(nodes_array)


    fig = plt.gcf()
    ax = fig.gca()
    plt.scatter(nodes_array[:,2], nodes_array[:,1], c='r', marker='o')
    ax.set_title('All nodes on this cross-section')
    ax.set_ylabel('y (m)')
    ax.set_xlabel('z (m)')
    plt.show()

#nodes_cross_section(X_nodes[0])
#nodes_cross_section()
