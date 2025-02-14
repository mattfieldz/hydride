import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import time


length = 10
time_length = 3
interface = 1
nodes = 200

D_oxide = 1
D_metal = 10


dx = length/nodes
dy = length/nodes
dt =  0.5 * 0.5 * dx ** 2 / max([D_oxide,D_metal])

r_o = D_oxide * dt/(dx**2)
r_m = D_metal * dt/(dx**2)

interface_grid = int(nodes*interface/length) 

concentration = np.zeros((nodes,nodes))

concentration[:,0] = 1 + 10*np.exp(-np.linspace(-5,5,nodes)**2)
concentration[:,-1] = 0



fig,axis = plt.subplots()
pcm = axis.pcolormesh(concentration,cmap=plt.cm.jet,vmin = 0, vmax = 1)
plt.colorbar(pcm,ax=axis)

def explicit_differences(concentration):
    counter = 0 
    while counter < time_length:
        copy_conc = np.copy(concentration)
        for i in range(1,nodes-1):
            
            for j in range(1,interface_grid):
            
                ddx = (copy_conc[i-1,j] - 2*copy_conc[i,j] + copy_conc[i+1,j])/(dx**2)
                ddy = (copy_conc[i,j-1] - 2*copy_conc[i,j] + copy_conc[i,j+1])/(dy**2)

                concentration[i,j] = D_oxide * dt * (ddx + ddy) + copy_conc[i,j]
            for j in range(interface_grid, nodes-1):
                
                ddx = (copy_conc[i-1,j] - 2*copy_conc[i,j] + copy_conc[i+1,j])/(dx**2)
                ddy = (copy_conc[i,j-1] - 2*copy_conc[i,j] + copy_conc[i,j+1])/(dy**2)

                concentration[i,j] = D_metal * dt * (ddx + ddy) + copy_conc[i,j]


        counter += dt

        pcm.set_array(concentration)
        axis.set_title(counter)
        plt.pause(0.01)
    return concentration
def implicit_differences(concentration):
    counter = 0 
    n = nodes

    conc_array = np.zeros((n,n))

    conc_array[0,0] = 1
    conc_array[n-1,n-1] = 1
    conc_array[0,1] = 0 
    conc_array[1,0] = 0
    conc_array[n-2,n-1] = 0 
    conc_array[n-1,n-2] = 0

    
    for i in range(1,interface_grid):
        conc_array[i,i] = (1 + 2 * r_o) 
        conc_array[i,i-1] = -r_o 
        conc_array[i,i+1] = -r_o 

    for i in range(interface_grid, n-1):
        conc_array[i,i] = (1 + 2 * r_m) 
        conc_array[i,i-1] = -r_m 
        conc_array[i,i+1] = -r_m 
    # print(conc_array)
    while counter < time_length:
        copy_conc = np.copy(concentration)
        
       
        
        concentration = np.linalg.solve(conc_array,copy_conc)
        counter += dt
        # print(concentration)
        # pcm.set_array([concentration])
        # axis.set_title(counter)
        # plt.pause(0.01)

    return concentration

def implicit_differences_sparse(concentration):
    counter = 0 
    n = nodes

    conc_array = np.zeros((3,nodes))
 
    

    for i in range(1,interface_grid):
        conc_array[0,i] = -r_o
        conc_array[2,i] = -r_o
        conc_array[1,i] = 1 + 2 * r_o

    for i in range(interface_grid,n-1):
        conc_array[0,i] = -r_m
        conc_array[2,i] = -r_m
        conc_array[1,i] = 1 + 2 * r_m
    
    conc_array[1,0] = 1
    conc_array[1,n-1] = 1
    conc_array[0,1] = 0
    conc_array[0,n-1] = -r_o
    conc_array[2,0] = -r_m
    conc_array[2,n-2] = 0

    while counter < time_length:
        copy_conc = np.copy(concentration)
        
        concentration = sp.linalg.solve_banded((1,1),conc_array,copy_conc)
        concentration[0] = 1
        concentration[-1] = 0
        counter += dt
        # print(counter)
        # print(concentration)
        pcm.set_array([concentration])
        axis.set_title(counter)
        plt.pause(0.01)

    return concentration


if __name__ == '__main__':
    x = explicit_differences(concentration)