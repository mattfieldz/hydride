import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import time


length = 10
time_length = 3
interface = 5
nodes = 200

D_oxide = 1
D_metal = 1


dx = length/nodes
dt =  0.5 * dx ** 2 / max([D_oxide,D_metal])

r_o = D_oxide * dt/(dx**2)
r_m = D_metal * dt/(dx**2)

interface_grid = int(nodes*interface/length) 

concentration = np.zeros(nodes)

concentration[0] = 1



fig,axis = plt.subplots()
pcm = axis.pcolormesh([concentration],cmap=plt.cm.jet,vmin = 0, vmax = 1)
plt.colorbar(pcm,ax=axis)

def explicit_differences(concentration):
    counter = 0 
    while counter < time_length:
        copy_conc = np.copy(concentration)
        
        for n in range(1,interface_grid):
            concentration[n] = D_oxide * dt * (copy_conc[n-1] - 2*copy_conc[n] + copy_conc[n+1])/(dx**2) + copy_conc[n]
        for n in range(interface_grid, nodes-1):
            concentration[n] = D_metal * dt * (copy_conc[n-1] - 2*copy_conc[n] + copy_conc[n+1])/(dx**2) + copy_conc[n]


        counter += dt

        # pcm.set_array([concentration])
        # axis.set_title(counter)
        # plt.pause(0.01)
    return concentration

def explicit_difference_matrix(concentration):
    counter = 0
    n = nodes
    # conc_array = np.zeros((3,nodes))
    # conc_array[0,1:n] = D_oxide * dt / dx**2 * np.ones(n-1)
    # conc_array[1,1:n-1] = -2 * D_oxide * dt / dx**2 * np.ones(n-2)
    # conc_array[2,0:n-1] = D_oxide * dt / dx**2 * np.ones(n-1)
    # conc_array[1,0] = 1
    # conc_array[1,n-1] = 1
    a = np.zeros(n)
    a[0] = 1

    

    while counter < time_length:
        copy_conc = np.copy(concentration)
        concentration = conc_array @ copy_conc + copy_conc
        print(concentration)

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

def crank_nicholson(concentration):
    counter = 0
    n = nodes
    
    data = np.array([])

    offsets = np.array([-1, 0, 1])

    a = sp.sparse.dia_array((data, offsets), shape=(n*n, n*n))

if __name__ == '__main__':
    # x = explicit_differences(concentration)
    # t0 = time.time()
    # y = implicit_differences(concentration)
    # t1 = time.time()
    # # print(x-y)
    # z = implicit_differences_sparse(concentration)
    # t2= time.time()
    # print(z)
    # print(y)
    # print(z-y)
    # print(t1-t0)
    # print(t2-t1)
    w = crank_nicholson(concentration)
    # explicit_difference_matrix(concentration)
