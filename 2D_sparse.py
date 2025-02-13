from matplotlib import pyplot as plt
import numpy as np
import scipy as sp

def matrix_generator(n,r):
    M=n**2
    diag = np.zeros(M)
    upper_diag = np.zeros(M)
    lower_diag = np.zeros(M)
    upper_diag_n = np.zeros(M)
    lower_diag_n = np.zeros(M)

    upper_diag1 = np.zeros(M)
    lower_diag1 = np.zeros(M)
    upper_diag_n1 = np.zeros(M)
    lower_diag_n1 = np.zeros(M)
    


    for i in range(M):

        # sets dirichlet surface condition
        if i < n:
            diag[i] = 1
        # sets neumann no flux conditions to second order accuracy
        elif i > M-n-1:
            diag[i] = -3
            lower_diag_n[i-n] = 4
            lower_diag_n1[i-n-1] = -1
        elif np.mod(i,n) == n-1:
            diag[i] = -3
            lower_diag[i-1] = 4
            lower_diag1[i-2] = -1
        elif np.mod(i,n) == 0:
            diag[i] = -3
            upper_diag[i+1] = 4
            upper_diag1[i+2] = -1
        # sets the rest of matrix with finite difference scheme
        else:
            diag[i] = 1 + 2*r[i]
            upper_diag[i+1] = -r[i]/2
            upper_diag_n[i+n] = -r[i]/2
            lower_diag[i-1] = -r[i]/2
            lower_diag_n[i-n] = -r[i]/2
    
    # Builds matrix as a sparse array
    data = np.array([lower_diag_n1, lower_diag_n ,lower_diag1,lower_diag  ,diag ,upper_diag,upper_diag1,upper_diag_n,upper_diag_n1 ])

    offsets = np.array([-n-1,-n,-2,-1, 0, 1,2,n,n+1])

    a = sp.sparse.dia_array((data, offsets), shape=(M, M))
    # Returns a 9x(N^2) array encoding all the non zero entries 
    return a.tocsr()

def c_vector_generator(n,r,concentration):
# Initialise

    c_matrix = np.zeros((n,n))
# Boundaries

    c_matrix[0,:] = 1 + 10*np.exp(-np.linspace(-5,5,n)**2)
    c_matrix[1:n,0] = 0
    c_matrix[1:n,n-1] = 0
    c_matrix[n-1,:] = 0
# Crank Nicholson terms
    for i in range(1,n-1):
            for j in range(1,n-1):
                c_matrix[i,j] = (1-2*r[i,j])*concentration[i,j] + 1/2 * r[i,j] * (concentration[i,j-1] + concentration[i,j+1] + concentration[i+1,j] + concentration[i-1,j])
    return c_matrix.flatten()



counter = 0 
total_time = 100
n=5

M = n*n

L = 10
interface = 1

dx = L/n
dt = 0.1

D_oxide = 0.01
D_metal = 1
D_hydride = 10**(-6)
diffusivity_matrix = np.zeros((n,n))


interface_grid_value = int(interface/L * n)

oxide = np.zeros((n,n))
metal = np.zeros((n,n))
hydride = np.zeros((n,n))

oxide[0:interface_grid_value,:] = 1
metal[interface_grid_value:n,:] = 1


concentration = np.zeros((n,n))

fig,axis = plt.subplots(2,2)
pcm_conc = axis[0,0].pcolormesh(concentration,cmap=plt.cm.jet,vmin = 0, vmax = 1)
# pcm_reaction = axis[1,1].pcolormesh(reaction,cmap=plt.cm.jet,vmin = 0, vmax = 1)
axis[0,0].set_title('H Concentration')
axis[1,0].set_title('Uranium Hydride')
axis[0,1].set_title('Uranium')
axis[1,1].set_title('Contour reaction term')
plt.colorbar(pcm_conc,ax=axis)

while counter < total_time:
    diffusivity_matrix = D_oxide * oxide + D_hydride * hydride + D_metal * metal
    r = diffusivity_matrix * dt / dx**2

    matrix = matrix_generator(n,r.flatten())

    plt.imshow(matrix.toarray().astype(bool))
    plt.show()

    c_vector = c_vector_generator(n,r,concentration)
    concentration_flatten = sp.sparse.linalg.spsolve(matrix,c_vector)
    concentration = concentration_flatten.reshape((n,n))
    
    counter += dt

    pcm_conc.set_array(concentration)
    # axis.set_title(counter)
    fig.suptitle(counter)
    plt.pause(0.01)