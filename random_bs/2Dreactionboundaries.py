import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import time


length = 10
time_length = 10
interface = 0.5
nodes = 5

k = 10
conc_crit = 0.1


D_oxide = 1
D_metal = 1
D_hydride = 0.01

dx = length/nodes
dy = length/nodes
dt = dx ** 2 / max([D_oxide,D_metal,D_hydride])
print('dt',dt)
r_o = D_oxide * dt/(dx**2)
r_m = D_metal * dt/(dx**2)

interface_grid = int(nodes*interface/length) 

concentration = np.zeros((nodes,nodes))
uranium = np.zeros((nodes,nodes))
uranium_hydride = np.zeros((nodes,nodes))
oxide = np.zeros((nodes,nodes))
reaction = np.zeros((nodes,nodes))
uranium[:,interface_grid:nodes] = 1
oxide[:,0:interface_grid] = 1

concentration[:,0] = 1 + 10*np.exp(-np.linspace(-5,5,nodes)**2)

def reac_func(conc, conc_crit,uran):
    reac = k * np.heaviside(conc-conc_crit,0.5) * conc**3 * uran
    return reac

x = y = np.linspace(-5,5,nodes)
X, Y = np.meshgrid(x,y)

fig,axis = plt.subplots(2,2)
pcm_conc = axis[0,0].pcolormesh(concentration,cmap=plt.cm.jet,vmin = 0, vmax = 1)
pcm_uranium = axis[0,1].pcolormesh(uranium,cmap=plt.cm.jet,vmin = 0, vmax = 1)
pcm_uranium_hydride = axis[1,0].pcolormesh(uranium_hydride,cmap=plt.cm.jet,vmin = 0, vmax = 1)
# pcm_reaction = axis[1,1].pcolormesh(reaction,cmap=plt.cm.jet,vmin = 0, vmax = 1)
axis[0,0].set_title('H Concentration')
axis[1,0].set_title('Uranium Hydride')
axis[0,1].set_title('Uranium')
axis[1,1].set_title('Contour reaction term')
plt.colorbar(pcm_conc,ax=axis)

def explicit_differences(concentration):
    counter = 0 
    while counter < time_length:
        copy_conc = np.copy(concentration)

        reaction = np.zeros((nodes,nodes))
        copy_conc[0,:] = concentration[1,:]
        copy_conc[-1,:] = concentration[-2,:]
        copy_conc[:,-1] = concentration[:,-2]

        for i in range(1,nodes-1):
            
            for j in range(1,nodes-1):
            
                D = D_metal * uranium[i,j] + D_oxide * oxide[i,j] + D_hydride * uranium_hydride[i,j]
                ddx = (copy_conc[i-1,j] - 2*copy_conc[i,j] + copy_conc[i+1,j])/(dx**2)
                ddy = (copy_conc[i,j-1] - 2*copy_conc[i,j] + copy_conc[i,j+1])/(dy**2)

                concentration[i,j] = dt * (D * (ddx + ddy) - 3 * reac_func(copy_conc[i,j],conc_crit,uranium[i,j])) + copy_conc[i,j]
                uranium_hydride[i,j] = dt * reac_func(copy_conc[i,j],conc_crit,uranium[i,j]) + uranium_hydride[i,j]
                uranium[i,j] = - dt * reac_func(copy_conc[i,j],conc_crit,uranium[i,j]) + uranium[i,j]

        counter += dt

        pcm_conc.set_array(concentration)
        pcm_uranium.set_array(uranium)
        pcm_uranium_hydride.set_array(uranium_hydride)
        axis[1,1].contour(X,Y,reac_func(concentration,conc_crit,uranium))
        fig.suptitle(counter)
        # axis.set_xlim((0,10))
        # axis.set_ylim((0,10))
        plt.pause(0.01)
    return concentration
def implicit_differences(concentration):
    conc_copy = np.copy(concentration)


    counter = 0 
    n = nodes
    m = nodes
    c_vector = np.zeros((n,m))

    

    conc_array = np.zeros((n*m,n*m))
    ex = np.ones(n*m)
    

    r = dt / (dx**2)

    ex = np.ones(n**2)

    exdiag = np.copy(ex)

    exdiag = (1 + 2*r) * exdiag


    ex1 = np.copy(ex)
    ex2 = np.copy(ex)
    ex3 = np.copy(ex)
    ex4 = np.copy(ex)

    for count in range(n):
        ex1[n*(count) - n] = 0
        ex2[n*(count)-1] = 0
        ex3[n*(count)+1] = 0
        ex4[n*(count) + - 2 * n] = 0
        ex3[n*(count)] = 0

        exdiag[n*(count)] = 1
        exdiag[n*(count)-1] += -r/2
        
        

    for count in range(n-1):
        exdiag[(count)+1] += -r/2 
        exdiag[-((count+1))] += -r/2
    


    data = np.array([-r * ex1, -r * ex2, exdiag, -r * ex3, -r * ex4])

    offsets = np.array([-n,-1, 0, 1,n])

    a = sp.sparse.dia_array((data, offsets), shape=(n*n, n*n))

    print(a.toarray())
    # print(a.toarra?y())
    # print(conc_array)
    while counter < time_length:
        copy_conc = np.copy(concentration)
        
        for i in range(0,n):
            for j in range(0,m):
                if j==0:
                    c_vector[i,j] = copy_conc[i,j]
                elif j==n-1 and i != n-1:
                    c_vector[i,j] = (1-2*r+r/2)*copy_conc[i,j] + r/2*(copy_conc[i-1,j] + copy_conc[i+1,j] + copy_conc[i,j-1])
                elif i==0 and j != n-1:
                    c_vector[i,j] = (1-2*r+r/2)*copy_conc[i,j] + r/2*(copy_conc[i+1,j] + copy_conc[i,j+1] + copy_conc[i,j-1])
                elif i==n-1 and j != n-1:
                    c_vector[i,j] = (1-2*r+r/2)*copy_conc[i,j] + r/2*(copy_conc[i-1,j] + copy_conc[i,j+1] + copy_conc[i,j-1])
                elif j == n-1 or i==n-1:
                    c_vector[i,j] = 0
                else:
                    c_vector[i,j] = (1-2*r)*copy_conc[i,j] + r/2*(copy_conc[i-1,j] + copy_conc[i+1,j] + copy_conc[i,j+1] + copy_conc[i,j-1])
        # print(c_vector)
        c = c_vector.flatten()
        # print(c)
        concentration_flattened = sp.sparse.linalg.spsolve(a.tocsr(),c)
        concentration = concentration_flattened.reshape((n,m))
        # print(concentration,'conc')
        # print(concentration_flattened)
        counter += dt
        # print(concentration)
        pcm_conc.set_array(concentration)
        # axis.set_title(counter)
        fig.suptitle(counter)
        plt.pause(0.01)
        
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

def implicit_not_sparse(concentration):
    conc_copy = np.copy(concentration)


    counter = 0 
    n = nodes
    m = nodes
    c_vector = np.zeros((n,m))

    M = n*n

    conc_array = np.zeros((n*m,n*m))
   
    r = dt / (dx**2)

   

    for i in range(M):
        for j in range(M):
            if i < n:
                conc_array[i,j] = 1
            elif i > M-n-1:
                conc_array[i,j] = 1
            elif np.mod(i,n):
    


   

    # print(a.toarra?y())
    # print(conc_array)
    while counter < time_length:
        copy_conc = np.copy(concentration)
        
        for i in range(0,n):
            for j in range(0,m):
                if j==0:
                    c_vector[i,j] = copy_conc[i,j]
                elif j==n-1 and i != n-1:
                    c_vector[i,j] = (1-2*r+r/2)*copy_conc[i,j] + r/2*(copy_conc[i-1,j] + copy_conc[i+1,j] + copy_conc[i,j-1])
                elif i==0 and j != n-1:
                    c_vector[i,j] = (1-2*r+r/2)*copy_conc[i,j] + r/2*(copy_conc[i+1,j] + copy_conc[i,j+1] + copy_conc[i,j-1])
                elif i==n-1 and j != n-1:
                    c_vector[i,j] = (1-2*r+r/2)*copy_conc[i,j] + r/2*(copy_conc[i-1,j] + copy_conc[i,j+1] + copy_conc[i,j-1])
                elif j == n-1 or i==n-1:
                    c_vector[i,j] = 0
                else:
                    c_vector[i,j] = (1-2*r)*copy_conc[i,j] + r/2*(copy_conc[i-1,j] + copy_conc[i+1,j] + copy_conc[i,j+1] + copy_conc[i,j-1])
        # print(c_vector)
        c = c_vector.flatten()
        # print(c)
        concentration_flattened = sp.sparse.linalg.spsolve(a.tocsr(),c)
        concentration = concentration_flattened.reshape((n,m))
        # print(concentration,'conc')
        # print(concentration_flattened)
        counter += dt
        # print(concentration)
        pcm_conc.set_array(concentration)
        # axis.set_title(counter)
        fig.suptitle(counter)
        plt.pause(0.01)
        


if __name__ == '__main__':
    # x = explicit_differences(concentration)
    # y = implicit_differences(concentration)
    z = implicit_not_sparse(concentration)