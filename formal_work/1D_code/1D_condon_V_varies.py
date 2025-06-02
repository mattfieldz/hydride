
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

n2 = 10001  # metal
# number of variables in the bulk
nv = 2  # [H] and [U] in the bulk plus diffusivity of mixture

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 1e3 * 1e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 20 * 1.0e-7  # oxide domain lengthscale 10nm
D1star = 1.0e-13  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 1.49e-10  # cm^2/s, diffusion of H in U (room temp value)
D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 (room temp value)
N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 1.0e-5  # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castar = 1.0e-4  # mol/cm^3, surface value for [H] -- *ad-hoc* value

# fixed reference values to keep time/lengthscales consistent in comparisons
Lrefstar = 100 * 1e-7  # using 1um=1000nm here
# Lrefstar = L2star

Drefstar = D2star  # using the U value

# non-dimensional domains
L3 = L3star / Lrefstar  # oxide
L2 = L2star / Lrefstar  # bulk


# non-dimensional max = 10^4 sec
tmax = 3.0e6 * Drefstar / (Lrefstar * Lrefstar)

# non-dimensional variables


D3 = D3star / Drefstar
D2 = D2star / Drefstar
D1 = D1star / Drefstar
cs = Csstar / Castar  # saturation value \in [0,1)
eps = Castar / N2star  # relative forcing measure

kstar = 1e14
reactK = (
    Lrefstar**2*eps**2*N2star**3/Drefstar * kstar  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
)
D = D2 * np.ones(n2, dtype=np.float64)  # mixture diffusion in the bulk, initially =D2

# information output sanity check
print(f"Diffusivity ratio D1*/Dref*={D1star/Drefstar}")
print(f"Diffusivity ratio D2*/Dref*={D2star/Drefstar}")
print(f"Diffusivity ratio D3*/Dref*={D3star/Drefstar}")
print(f"Lengthscale ratio L2*/Lref*={L2star/Lrefstar}")
print(f"Lengthscale ratio L3*/Lref*={L3star/Lrefstar}")
print(f"Diffusion timescale for L1* of oxide T3*={L3star**2/D3star} sec")
print(f"Diffusion timescale for Lref* of metal T2*={Lrefstar**2/D2star} sec")
print(f"Diffusion timescale for Lref* of hydride T1*={Lrefstar**2/D1star} sec")
print(f"Ratio of H to U concentration eps={eps}")
print(f"Relative solubility limit ={cs}")
print(f"Non-dimensional max time ={tmax}")
print(f"Non-dimensional reaction constant ={reactK}")
# metric_file = open("1Dexamples/data/metric.dat", "w")

# step sizes

h2 = 0.001
h2s = h2 * h2
dt = 0.0001
tol = 1.0e-8
alpha_c = 0.98



m = 1 # reaction order

# time step the system
t = 0
step = 1  # periodic output counter
iteration = 0
# Newton iteration loop


c2 = np.exp(-np.linspace(0,100,n2)) * np.ones(n2, dtype=np.float64)  # metal
alpha = alpha_c * np.ones(n2, dtype=np.float64)  # volume fraction of metal
c_xi_0_list = []

alpha_c_values = [0.9+0.005*i for i in range(20)]
m_values = [1+0.4*i for i in range(10)]
cs_values = [0.01+0.01*i for i in range(10)]
# cs_values = [0.000116]


D2star = 2.7652 * 1e-10
kstar = 15481


length_scale = (D2star / kstar) ** 0.5

domain = n2 * h2
c_xi_computed_values = []
print(length_scale,domain)

# rough_estimate = [1.7, 1.7, 1.7, 1.7, 1.6, 1.6, 1.6, 1.6, 1.6, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.4, 1.4, 1.4, 1.4, 1.4]
# rough_estimate = [1.7, 1.69, 1.67, 1.65, 1.63, 1.62, 1.6, 1.58, 1.57, 1.54, 1.53, 1.51, 1.5, 1.48, 1.46, 1.44, 1.43, 1.41, 1.39, 1.38]
# rough_estimate = [1.703, 1.686, 1.669, 1.652, 1.634, 1.617, 1.6, 1.583, 1.566, 1.548, 1.531, 1.514, 1.497, 1.48, 1.462, 1.445, 1.428, 1.411, 1.394, 1.376]
rough_estimate = [1.7032, 1.686, 1.6688, 1.6516, 1.6344, 1.6172, 1.6, 1.5828, 1.5656, 1.5484, 1.5312, 1.514, 1.4968, 1.4796, 1.4624, 1.4452, 1.428, 1.4108, 1.3936, 1.3764]

# print('lol',np.array(rough_estimate)/(1-np.array(cs_values)))

m_val_array = np.array(m_values)


rough_estimate = np.sqrt(6/(m_val_array+1))
print(rough_estimate)
rough_estimate = [1.7204, 1.5712, 1.4564, 1.3656, 1.2914, 1.2292, 1.176, 1.1298, 1.0892, 1.0532]

print(rough_estimate)

estimate_counter = 0
for m in m_values:
    
    c_xi_values = [rough_estimate[estimate_counter]-0.00008+0.00002*i for i in range(10)]

    estimate_counter +=1

    min_value_list = []
    min_value_dictionary = []
    for i in (c_xi_values):
        V = i/(3*(1-alpha_c))
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
                    b += [1-c2[j],alpha_c-alpha[j]]
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
                    
                    row += 4 * [k+1]
                    col += [k-nv+1,k,k+1,k+nv+1]
                    val += [
                        1/(2*h2) * V,
                        m * c2[j]**(m-1) * alpha[j],
                        c2[j]**m,
                        -1/(2*h2) * V,
                    ]
                    b += [a_xi * V - c2[j]**m * alpha[j]]
            
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

            
            if iteration > 10:
                print(residual)
                
                # print("Too many iterations")
                # sys.exit()
            # print(residual)
        # c_xi_0_list.append(c_xi_0)
        # print(cs)
    # print(c_xi_0_list)
        min_value_dictionary.append(i)
        min_value_list.append(-c_xi_0/(3*(1-alpha_c))-V)
        print(i,-c_xi_0/(3*(1-alpha_c))-V)
    min_value_array = np.array(min_value_list)
    
    c_xi_computed_values.append(round(min_value_dictionary[np.argmin(np.abs(min_value_array))],5))
    print(c_xi_computed_values)

# np.savetxt(
#     f"formal_work/data1D/c_xi_0_list_cs_computed_varying.dat",
#     np.transpose(np.array([cs_values, -np.array(c_xi_computed_values)])),
# )

print(alpha)

length_scale = (D2star / kstar) ** 0.5

domain = n2 * h2

interpolated_alpha = sp.interpolate.CubicSpline(np.linspace(0,domain,n2)*length_scale, alpha - 0.999)
roots = interpolated_alpha.roots(extrapolate=False)
print('roots = ',roots)

# fig, axes = plt.subplots(2)
# axes[0].plot(np.linspace(0,domain,n2)*length_scale * 1e7,alpha)
# axes[1].plot(np.linspace(0,domain,n2)*length_scale * 1e7,c2+cs)

# axes[0].set_ylabel(r'$\alpha$')
# axes[0].set_xlabel(r'$x^* / \mathrm{nm}$')

# axes[1].set_ylabel(r'$c$')
# axes[1].set_xlabel(r'$x^* / \mathrm{nm}$')

plt.plot(m_values,rough_estimate,label='Numerical solution')
plt.plot(m_values,np.sqrt(6/(m_val_array+1)),label='Asymptotic solution')
plt.xlabel('Reaction order m')
plt.ylabel(r'$C_\xi(0)$')

plt.legend()

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


