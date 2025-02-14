
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
# import deps.sp2petsc as petsc

# indexing is 1=Hydride (not present), 2=Metal, 3=Oxide
# discretisation
n3 = 101  # oxide points
n2 = 1001  # metal points

# e1 = 4 to jump to concentration nodes
e1 = 4

# dimensional quantities, natchiar 2024 for diffusivity of hydrogen in U, UH3, UO2
L3star = 1.0e-6  # oxide lengthscale 10 nm = 1.e-8 m = 1.e-6 cm
L2star = 5.0e-3  # bulk lengthscale  50 um = 5.e-5 m = 5.e-3 cm
D3star = 1.18e-13  # cm^2/s, diffusion of H in UO2 @ room temp
D2star = 1.49e-10  # cm^2/s, diffusion of H in U @ room temp
D1star = 8.94e-16


# non-dimensional domains
L3 = L3star / L2star
L2 = 1
x3_nodes = np.linspace(0, L3, n3)
x2_nodes = np.linspace(L3, L3 + L2, n2)
tmax = 10

# information sanity check
print(f"Diffusivity ratio D3*/D2*={D3star/D2star}")
print(f"Lengthscale ratio L3*/L2*={L3star/L2star}")
print(f"Diffusion time scale for oxide T1*={L3star**2/D3star} sec")
print(f"Diffusion time scale for metal T2*={L2star**2/D2star} sec")

# non-dimensional variables
c3_nodes = np.zeros(n3*e1, dtype=np.float64)
c2_nodes = np.zeros(n2*e1, dtype=np.float64)
D3 = D3star / D2star
D1 = D1star / D2star
D2 = 1

# step sizes
h3 = L3 / (n3 - 1)
h3s = h3 * h3
h2 = L2 / (n2 - 1)
h2s = h2 * h2
dt = 0.001

# time step
t = 0
step = 1

# Plotting
fig,axis = plt.subplots(2,2)



while t < tmax:
    row = []
    col = []
    val = []
    b = []  # RHS vector for Ax=b

    # surface has C3=1 for j=0
    row += [0,1,2,3]
    col += [0,1,2,3]
    val += [1,1,1,1]
    b += [1,0,0,1]

    # surface a1 = 0, a2 = 0, a3 = 1
    





    # internal nodes
    for j in range(1, n3 - 1):
        # concentration nodes
        row += 3 * [e1*j] # this is [j,j,j]
        col += [e1*j - e1, e1*j, e1*j + e1]
        val += [D3 / h3s, -2 * D3 / h3s - 1 / dt, D3 / h3s]
        b += [-c3_nodes[j] / dt]
        # metal nodes
        row += [e1 * j + 1, e1 * j + 2, e1 * j + 3]
        col += [e1 * j + 1, e1 * j + 2, e1 * j + 3]
        val += [1,1,1]
        b += [0,0,1]

    # interface has equal concentrations: C3 = C2
    j = n3 - 1
    row += 2 * [e1*j]
    col += [e1*j, e1*j + e1]
    val += [1, -1]
    b += [0]
# oxide region
    row += [e1 * j + 1, e1 * j + 2, e1 * j + 3]
    col += [e1 * j + 1, e1 * j + 2, e1 * j + 3]
    val += [1,1,1]
    b += [0,0,1]

    # interface has matched flux: D3*C3' - D2*C2' = 0 
    j = n3
    D = D2 * c2_nodes[2] + D1 * c2_nodes[1]
    row += 4 * [e1*j]
    col += [e1*(j - 2), e1*(j - 1), e1*j, e1*(j+1)]  # C3_{n3-2},C3_{n3-1},C2_{0},C2_{1}
    # note: only first order accurate, but can use 3 points either side for second order
    val += [-D3 / h3, D3 / h3, D / h2, -D / h2]
    b += [0]

    row += [e1 * j + 1, e1 * j + 2, e1 * j + 3]
    col += [e1 * j + 1, e1 * j + 2, e1 * j + 3]
    val += [1,1,1]
    b += [0,1,0]

    # internal nodes in metal region
    for j in range(1, n2 - 1):  # logical node in metal
        k = j + n3  # offset by oxide region eqns

        D = D2 * c2_nodes[2+e1*j] + D1 * c2_nodes[1+e1*j]

        row += 3 * [e1 * k]
        col += [e1 * k - e1, e1*k, e1*k + e1]
        val += [D / h2s, -2 * D / h2s - 1 / dt, D / h2s]
        b += [-c2_nodes[e1 * j] / dt]

        row += [e1 * k + 1, e1 * k + 2, e1 * k + 3]
        col += [e1 * k + 1, e1 * k + 2, e1 * k + 3]
        val += [1,1,1]
        b += [0,1,0]

    # interior boundary has zero first order deriv 
    j = n2 - 1  # logical node in metal
    k = j + n3  # row in mtx
    row += 2 * [e1 * k]
    col += [e1 * k, e1 * k - e1]
    val += [1 / h2, -1 / h2]
    b += [0]

    row += [e1 * k + 1, e1 * k + 2, e1 * k+ 3]
    col += [e1 * k + 1, e1 * k + 2, e1 * k + 3]
    val += [1,1,1]
    b += [0,1,0]

    a = sp.sparse.coo_matrix((val, (row, col)), shape=(e1 * (n3 + n2), e1 * (n3 + n2)))

    # system = petsc.PETScSparseLinearSystem(a, b)

    # x = system.linear_solve()
    
    x = sp.sparse.linalg.spsolve(a.tocsr(),b)

    # store the solution for this time level
    c3_nodes[0:n3*e1] = x[0:n3*e1]
    c2_nodes[0:n2*e1] = x[n3*e1 : n3*e1 + n2*e1]
    t += dt
    # save every 10 steps
    if step % 10 == 0:
        print(
            f"Save at t={t} c_int={c3_nodes[e1*(n3-1)]}={c2_nodes[0]} c_end={c2_nodes[n2-1]}"
        )
        # np.savetxt(
        #     f"./c3_{t:.2f}.dat",
        #     np.transpose(np.array([x3_nodes, c3_nodes])),
        # )
        # np.savetxt(
        #     f"./c2_{t:.2f}.dat",
        #     np.transpose(np.array([x2_nodes, c2_nodes])),
        # )
        
    pcm_conc = axis[0,0].plot(x[::4])

    # axis.set_title(counter)
    fig.suptitle(step)
    plt.pause(0.01)
    step += 1
