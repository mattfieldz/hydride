import numpy as np
import scipy as sp
n = m = 4
M = n**2
ex = np.zeros(M)
dx = 1
dt = 2

r = dt / (dx**2)

ex2l = np.copy(ex)
ex1l = np.copy(ex)
ex1u = np.copy(ex)
ex2u = np.copy(ex)

for i in range(M):
    if i < n:
        ex[i] = 1
    elif i > M-n-1:
        ex[i] = 1
    else:
        ex[i] = 2
        





data = np.array([ex2l ,ex1l  ,ex ,ex1u,ex2u ])

offsets = np.array([-n,-1, 0, 1,n])

a = sp.sparse.dia_array((data, offsets), shape=(n*m, n*m))

print(a.toarray())