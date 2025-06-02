
"""e.g. python diffusionBC_reaction.py -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist"""

import numpy as np
import scipy as sp
# import deps.sp2petsc as petsc
import sys
from matplotlib import pyplot as plt

# some functions to remap coords
def x_from_X(X):
    """Defines the (non-uniformly distributed) physical coordinate 'x'
    in terms of the (uniformly distributed) computational coordinate
    'X'. Returns: not only x(X) but also x'(X) and x''(X).

    """
    BX = 0.01  # smaller BX leads to more NON-uniformity
    # +0*X below prevents numpy from returning a scalar when acting
    # on an iterable container for X
    return (X + BX) * (X + BX) - BX * BX, 2 * (X + BX)


# indexing is 1=Hydride (not present), 2=Metal, 3=Oxide
#
# discretisation

n2 = 8001  # metal
# number of variables in the bulk
nv = 2  # [H] and [U] in the bulk plus diffusivity of mixture

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 1e2 * 1e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 20 * 1.0e-7  # oxide domain lengthscale 10nm
D1star = 1.0e-13  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 1.49e-10  # cm^2/s, diffusion of H in U (room temp value)
D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 (room temp value)
N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 1.0e-5  # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castar = 1.0e-4  # mol/cm^3, surface value for [H] -- *ad-hoc* value

# fixed reference values to keep time/lengthscales consistent in comparisons
# Lrefstar = 100 * 1e-7  # using 1um=1000nm here
# Lrefstar = L2star

reactK = Castar * N2star

Lrefstar = D2star**(0.5)/reactK**0.5


Drefstar = D2star  # using the U value


D3 = D3star / Drefstar
D2 = D2star / Drefstar
D1 = D1star / Drefstar
cs = Csstar / Castar  # saturation value \in [0,1)
eps = Castar / N2star  # relative forcing measure

# kstar = 1e14
# reactK = (
#     Lrefstar**2*eps**2*N2star**3/Drefstar * kstar  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
# )
D = D2 * np.ones(n2, dtype=np.float64)  # mixture diffusion in the bulk, initially =D2


# step sizes

length = 5


h2 = length/n2
h2s = h2 * h2
# dt = 0.0001
tol = 1.0e-8
alpha_c = 0.98

m = 1 # reaction order

# time step the system
t = 0
step = 1  # periodic output counter
iteration = 0
# Newton iteration loop


# c2 = np.exp(-np.linspace(0,100,n2)) * np.ones(n2, dtype=np.float64)  # metal
# alpha = (1-np.linspace(0,1,n2)) * alpha_c * np.ones(n2,dtype=np.float64) + np.linspace(0,1,n2) * np.ones(n2, dtype=np.float64)  # volume fraction of metal
c_xi_0_list = []

alpha_c_values = [0.9+0.005*i for i in range(20)]
m_values = [1+i for i in range(0,3)]
# cs_values = [0.01+0.01*i for i in range(20)]
# cs_values = [0.1]

for m in m_values:
    residual = 1
    iteration = 0
    c2 = np.ones(n2, dtype=np.float64)  # metal
    alpha = np.ones(n2, dtype=np.float64)  # volume fraction of metal       
    while residual > tol:
        iteration += 1
        # construct the sparse matrix
        row = []
        col = []
        val = []
        b = []

        for j in range(n2):
            k = j * nv
            if j == 0:
                
                row += [k,k+1]
                col += [k,k+1]
                val += [1,1]
                b += [1-cs-c2[j],alpha_c-alpha[j]]
            elif j == n2-1:
                
                row += [k,k+1]
                col += [k,k+1]
                val += [1,1]
                b += [-c2[j],1-alpha[j]]

            else:
                
                # c equation
                row += 4 * [k]
                col += [k-nv,
                        k,
                        k+nv,
                        k+1]
                val += [1/h2s,
                        -2/h2s - 3*m*c2[j]**(m-1) * alpha[j],
                        1/h2s,
                        -3*c2[j]**m
                ]
                b += [-(c2[j-1]-2*c2[j]+c2[j+1])/h2s + 3*c2[j]**m*alpha[j]]
                
                # alpha equation

                c_xi_0 = (-3*c2[0]+4*c2[1]-c2[2])/(2*h2)
                a_xi = (alpha[j+1]-alpha[j-1])/(2*h2)
                if j == 1:
                    row += 6 * [k+1]
                    col += [k-nv,
                            k-nv+1,
                            k,
                            k+1,
                            k+nv,
                            k+nv+1]
                    val += [
                        -3/(2*h2) * a_xi / (3*(1-alpha_c)),
                        -1/(2*h2) * c_xi_0 / (3*(1-alpha_c)),
                        4/(2*h2) * a_xi / (3*(1-alpha_c)) + m * c2[j]**(m-1) * alpha[j],
                        c2[j]**m,
                        -1/(2*h2) * a_xi / (3*(1-alpha_c)),
                        1/(2*h2) * c_xi_0 / (3*(1-alpha_c))
                    ]
                    b += [-a_xi * c_xi_0 /(3*(1-alpha_c)) - c2[j]**m * alpha[j]]
                elif j == 2:
                    row += 6 * [k+1]
                    col += [0,k-nv,k-nv+1,k,k+1,k+nv+1]
                    val += [
                        -3/(2*h2) * a_xi / (3*(1-alpha_c)),
                        4/(2*h2) * a_xi / (3*(1-alpha_c)),
                        -1/(2*h2) * c_xi_0 / (3*(1-alpha_c)),
                        -1/(2*h2) * a_xi / (3*(1-alpha_c)) + m * c2[j]**(m-1) * alpha[j],
                        c2[j]**m,
                        1/(2*h2) * c_xi_0 / (3*(1-alpha_c)),
                        
                    ]
                    b += [-a_xi * c_xi_0 /(3*(1-alpha_c)) - c2[j]**m * alpha[j]]
                else:
                    row += 7 * [k+1]
                    col += [0,2,4,k-nv+1,k,k+1,k+nv+1]
                    val += [
                        -3/(2*h2) * a_xi / (3*(1-alpha_c)),
                        4/(2*h2) * a_xi / (3*(1-alpha_c)),
                        -1/(2*h2) * a_xi / (3*(1-alpha_c)),
                        -1/(2*h2) * c_xi_0 / (3*(1-alpha_c)),
                        m * c2[j]**(m-1) * alpha[j],
                        c2[j]**m,
                        1/(2*h2) * c_xi_0 / (3*(1-alpha_c)),
                    ]
                    b += [-a_xi * c_xi_0 /(3*(1-alpha_c)) - c2[j]**m * alpha[j]]
        
        a = sp.sparse.coo_matrix((val, (row, col)), shape=(n2 * nv,n2 * nv))
        # system = petsc.PETScSparseLinearSystem(a, b)
        # x = system.linear_solve()

        # print('det',np.linalg.det(a.toarray()))
        
        x = sp.sparse.linalg.spsolve(a.tocsr(),b)
        # residual is the maximum correction value
        residual = sp.linalg.norm(x, ord=np.inf)
       
        # print(
        #    f"iteration = {iteration} residual = {residual} ||b||_inf = {sp.linalg.norm(b, ord=np.inf)}"
        # )
        # add the corrections to the current guess
        
        c2[0:n2] += x[0 : n2 * nv : nv]
        alpha[0:n2] += x[1 : n2 * nv : nv]

        
        # if iteration > 10:
            # print("Too many iterations")
            # sys.exit()
        print(residual)
    c_xi_0_list.append(c_xi_0)
    
print(c_xi_0_list)



# np.savetxt(
#     f"formal_work/data1D/c_xi_0_list_m_new_varying.dat",
#     np.transpose(np.array([m_values, c_xi_0_list])),
# )

print(alpha)

# plt.plot(m_values,c_xi_0_list)
# plt.plot(m_values,-np.sqrt(6/(np.array(m_values)+1)))

plt.plot(np.linspace(0,length,n2),c2)

plt.show()



#
# Solution for this time step has converged
#

# measure the interface value of [H] and position of [alpha]=0.01
# if step % 10 == 0:
#     hydride_interface = 0
    # find positions where [UH3] concentration is 99% in the bulk
    # interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - threshold_alpha)
    # roots = interpolated_alpha.roots(extrapolate=False)
    # if len(roots) > 0:
    #     hydride_interface = roots[0]
    # OUTPUT:
    # non-dim time, dim time (sec), non-dim [H] @interface, non-dim [U] @interface, dim UH3 depth (nm)
    # metric_file.write(
    #     f"{t} {t*(Lrefstar**2)/(Drefstar)} {c2[0]} {alpha[0]} {hydride_interface*Lrefstar*1.e7}\n"
    # )
    # metric_file.flush()


# np.savetxt(
#     f"formal_work/data1D/k1e17/c2{t:.2f}.dat",
#     np.transpose(np.array([x2_nodes, c2])),
# )


