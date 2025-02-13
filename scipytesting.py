import numpy as np
import scipy as sp

from scipy.sparse import dia_array
from matplotlib import pyplot as plt

n = 10
r = 2

ex = np.ones(n**2)
exdiag = np.copy(ex)

exdiag = (1 + 2*r) * exdiag


ex1 = np.copy(ex)
ex2 = np.copy(ex)
ex3 = np.copy(ex)
ex4 = np.copy(ex)

for count in range(n):
    # ex1[n*(count) - n] = 0
    ex2[n*(count)-1] = 0
    ex3[n*(count)+1] = 0
    # ex4[n*(count) + - 2 * n] = 0
    ex3[n*(count)] = 0

    # exdiag[n*(count)] = 1
#     exdiag[n*(count)-1] += -r/2
    
    

# for count in range(n-1):
#     exdiag[(count)+1] += -r/2 
#     exdiag[-((count+1))] += -r/2



data = np.array([-r * ex1, -r * ex2, exdiag, -r * ex3, -r * ex4])

offsets = np.array([-n,-1, 0, 1,n])

a = dia_array((data, offsets), shape=(n*n, n*n))

b = np.linspace(1,9,9)

b_mat = b.reshape((3,3))
print(b_mat)

print(b.flatten())

# print(a@b)
print(a.toarray())
# d = a.tocsr()
# c = sp.sparse.linalg.spsolve(d,b)


# print(a@c)
# print(a.toarray())

print(a.toarray()[0,1])

counter = 0 
total_time = 10
n=4
m = n
c_vector = np.zeros((n,m))

M = n*n

conc_array = np.zeros((n*m,n*m))


L = 10
dx = L/n
dt = 0.1

r = dt/dx**2

concentration = np.zeros((n,n))
r= 4
for i in range(M):
    for j in range(M):
        if i < n:
            conc_array[i,i] = 1
        elif i > M-n-1:
            conc_array[i,i] = 1
            conc_array[i,i-n] = -1
        elif np.mod(i,n) == n-1:
            conc_array[i,i] = 1
            conc_array[i,i-1] = -1
        elif np.mod(i,n) == 0:
            conc_array[i,i] = 1
            conc_array[i,i+1] = -1
        else:
            conc_array[i,i] = 1+2*r
            conc_array[i,i+1] = -r/2
            conc_array[i,i-1] = -r/2
            conc_array[i,i+n] = -r/2
            conc_array[i,i-n] = -r/2

c_matrix = np.zeros((n,n))


c_matrix[0,:] = 1 + 10*np.exp(-np.linspace(-5,5,n)**2)
c_matrix[1:n,0] = 0
c_matrix[1:n,n-1] = 0
c_matrix[n-1,:] = 0


print(conc_array,'conc')


fig,axis = plt.subplots(2,2)
pcm_conc = axis[0,0].pcolormesh(concentration,cmap=plt.cm.jet,vmin = 0, vmax = 1)
# pcm_reaction = axis[1,1].pcolormesh(reaction,cmap=plt.cm.jet,vmin = 0, vmax = 1)
axis[0,0].set_title('H Concentration')
axis[1,0].set_title('Uranium Hydride')
axis[0,1].set_title('Uranium')
axis[1,1].set_title('Contour reaction term')
plt.colorbar(pcm_conc,ax=axis)

while counter < total_time:

    for i in range(1,n-1):
        for j in range(1,n-1):
            c_matrix[i,j] = (1-2*r)*concentration[i,j] + r/2 * (concentration[i,j-1] + concentration[i,j+1] + concentration[i+1,j] + concentration[i-1,j])
    

    concentration_flatten = np.linalg.solve(conc_array,c_matrix.flatten())
    concentration = concentration_flatten.reshape((n,n))
    
    pcm_conc.set_array(concentration)
    # axis.set_title(counter)
    fig.suptitle(counter)
    plt.pause(0.01)



# print(conc_array)
