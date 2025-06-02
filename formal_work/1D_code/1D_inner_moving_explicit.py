
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


c2_tot = np.loadtxt(
        f'formal_work/data1D/c2_condon_oxide_larger_D{13.83:.2f}.dat',
        )
alpha_tot = np.loadtxt(
    f'formal_work/data1D/alpha_condon_oxide_larger_D{13.83:.2f}.dat',
    )




n2 = 101  # metal
# number of variables in the bulk
nv = 3  # [H] and [U] in the bulk plus diffusivity of mixture

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
Lrefstar = L2star

D_factor = 1e3
D_large = 1
# Drefstar = D_factor * D2star # using the U value
Drefstar = D2star * D_factor
# non-dimensional domains
L3 = L3star / Lrefstar  # oxide
L2 = L2star / Lrefstar  # bulk
# Choose computation coord Xmax such that x_from_X(Xmax)=L2
Xmax = sp.optimize.root(lambda X: x_from_X(X)[0] - L2, 5).x[0]
# exit if sanity check fails
assert abs(x_from_X(Xmax)[0] - L2) < 1.0e-8
# nodes in the oxide and bulk

X2_nodes = np.linspace(0, Xmax, n2)  # bulk nodes in comp. coordinate
# physical node locations
x2_nodes, x2d_nodes = x_from_X(X2_nodes)
# non-dimensional max = 10^4 sec
tmax = 3.0e6 * Drefstar / (Lrefstar * Lrefstar)

# Xmax = 1
# x2_nodes = np.linspace(0,100,n2)
# x2d_nodes = np.ones(n2)


# non-dimensional variables

c2 = 0.01*np.ones(n2, dtype=np.float64)  # metal
alpha = np.ones(n2, dtype=np.float64)  # volume fraction of metal
Gamma = 0.0


# c2[0] = 0.0001

D3 = D3star / Drefstar
D2 = D2star / Drefstar
D1 = D1star / Drefstar
# D1 = D2
cs = Csstar / Castar  # saturation value \in [0,1)
eps = Castar / N2star  # relative forcing measure






kstar = 1e4
m = 1   # reaction order

reactK = (
    Lrefstar**2*N2star/Drefstar * kstar  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
)



# tscaling = (1/(eps*reactK))
# Lrefstar = (D2star/reactK)**0.5

tscaling = Lrefstar**2/Drefstar


print('x scaling : ', Lrefstar)
print('t scaling : ', tscaling)


c2_spline = sp.interpolate.CubicSpline(c2_tot[:,0], (c2_tot[:,1]-cs ) * reactK**0.5)
alpha_spline = sp.interpolate.CubicSpline(alpha_tot[:,0],alpha_tot[:,1])

c2_gamma = sp.interpolate.CubicSpline(c2_tot[:,0], c2_tot[:,1]-cs)

roots = c2_gamma.roots(extrapolate=False)
if len(roots) > 0:
    Gamma_init = roots[0]

Gamma = reactK**0.5 * Gamma_init

x_nodes = np.linspace(0,Gamma_init,n2)

# c2 = c2_spline(x_nodes)
# alpha = alpha_spline(x_nodes)



print('roots',Gamma)


D = D2 * np.ones(n2, dtype=np.float64)  # mixture diffusion in the bulk, initially =D2

# information output sanity check
print(f"Diffusivity ratio Dref*={Drefstar}")
print(f"Lengthscale ratio L2*/Lref*={L2star/Lrefstar}")
print(f"Diffusion timescale for Lref* of hydride T1*={Lrefstar**2/D1star} sec")
print(f"Ratio of H to U concentration eps={eps}")
print(f"Relative solubility limit ={cs}")
print(f"Non-dimensional max time ={tmax}")
print(f"Non-dimensional reaction constant ={reactK}")
# metric_file = open("1Dexamples/data/metric.dat", "w")

# step sizes

h2 = 1/(n2-1)
h2s = h2 * h2
dt_standard = 0.00001

dt = dt_standard
tol = 1.0e-8
threshold_alpha = 0.98


Y_nodes = np.linspace(0,1,n2)

# time step the system
t = 0
step = 1  # periodic output counter

hydride_interface = []

t_0 = np.zeros(n2)

Gamma = 0.0000001

D = D2 * np.ones(n2)



beta_vector = np.zeros(n2*2+1)
calphagamma = np.zeros(n2*2+1)
calphagamma[0:n2] = c2
calphagamma[n2:2*n2] = alpha
calphagamma[2*n2] = Gamma


Matrix = np.zeros((2*n2+1,2*n2+1))

Gamma_list = []
t_list =[]
while t < tmax:
    step += 1
    t += dt
    D = D2 * alpha + D_large * (1-alpha)
    for j in range(n2):
       
        if j == 0:
            alpha_deriv = (-3 * alpha[j] + 4 * alpha[j+1] - alpha[j+2])/(2*h2)
            Aa = Gamma
            Ba = -Y_nodes[j] * alpha_deriv
            Ca =  Gamma * alpha[j] - Gamma * Y_nodes[j] * alpha_deriv - dt * Gamma * c2[j] * alpha[j]
            Matrix[n2+j,n2+j] = Aa
            Matrix[n2+j,2*n2] = Ba
            beta_vector[n2+j] = Ca

            # c condition

            Matrix[j,j] = -D[j]*3 / (2*h2)
            Matrix[j,j+1] = D[j]*(4) / (2*h2)
            Matrix[j,j+2] = -D[j]*1 / (2*h2)
            Matrix[j,2*n2] = -D3 * (cs-1)/(L3)
            beta_vector[j] = 0

            

        elif j == n2-1:
            alpha_deriv = (3 * alpha[j] - 4 * alpha[j-1] + alpha[j-2])/(2*h2)
            c_deriv = (3 * c2[j] - 4 * c2[j-1] + c2[j-2])/(2*h2)
             # C
            Matrix[j,j] = 1
            beta_vector[j] = 0
             # Alpha
            Matrix[n2+j,n2+j] = 1
            beta_vector[n2+j] = 1

            # Gamma condition
            
            Matrix[2*n2,j] = D[j]*3 / (2*h2)
            Matrix[2*n2,j-1] = D[j]*(-4) / (2*h2)
            Matrix[2*n2,j-2] = D[j]*1 / (2*h2)
            Matrix[2*n2,2*n2] = cs * ((D2*eps*reactK**0.5)/(np.pi*t))**0.5
            beta_vector[2*n2] = 0

        else:
            alpha_deriv = (alpha[j+1]-alpha[j-1]) / (2 * h2)
            c_deriv = (c2[j+1]-c2[j-1]) / (2 * h2)
            # c eqn

            Ac = eps * Gamma
            Bc = -eps * Y_nodes[j] * c_deriv
            Cc = c2[j] * eps * Gamma - eps * Gamma * Y_nodes[j] * c_deriv
            + dt * (D[j] * (c2[j-1] - 2 * c2[j] + c2[j+1]) / (h2s) - 3 * Gamma**2 * alpha[j] * c2[j])

            Aa = Gamma
            Ba = -Y_nodes[j] * alpha_deriv
            Ca = Gamma * alpha[j] - Gamma * Y_nodes[j] * alpha_deriv - dt * Gamma * c2[j] * alpha[j]

            Matrix[j,j] = Ac
            Matrix[j,2*n2] = Bc
            Matrix[n2+j,n2+j] = Aa
            Matrix[n2+j,2*n2] = Ba

            beta_vector[j] = Cc
            beta_vector[n2+j] = Ca
    
    # plt.imshow(Matrix.astype(bool))
    # plt.show()
    # print(beta_vector)
    # print(Matrix)
    calphagamma = np.linalg.solve(Matrix,beta_vector)
    c2 = calphagamma[0:n2]
    alpha = calphagamma[n2:2*n2]
    Gamma = calphagamma[2*n2]
    
    # print(c2)
    # print(alpha)
    
    if step % 100 == 0:
        t_list.append(t)
        Gamma_list.append(Gamma)
        print(alpha[0],t/(eps*reactK**0.5) * Lrefstar**2/Drefstar)
        plt.plot(np.array(t_list)/(eps*reactK**0.5) * Lrefstar**2/Drefstar,1e7 * Lrefstar * np.array(Gamma_list)/reactK**0.5)
        plt.pause(0.01)


    # Solution for this time step has converged
    #

    # measure the interface value of [H] and position of [alpha]=0.01
    # if step % 10 == 0:
        
        # find positions where [UH3] concentration is 99% in the bulk
        # interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - threshold_alpha)
        # interpolated_c2 = sp.interpolate.CubicSpline(x2_nodes, c2)
        # roots = interpolated_alpha.roots(extrapolate=False)
        # if len(roots) > 0:
        #     hydride_interface.append(roots[0])
            
           
        # else:
        #     hydride_interface.append(0)
        
        # OUTPUT:
        # non-dim time, dim time (sec), non-dim [H] @interface, non-dim [U] @interface, dim UH3 depth (nm)
        # metric_file.write(
        #     f"{t} {t*(Lrefstar**2)/(Drefstar)} {c2[0]} {alpha[0]} {hydride_interface*Lrefstar*1.e7}\n"
        # )
        # metric_file.flush()

    # save every 100 steps
    # if step % 100 == 0:
    #     # print(
    #     #     f"Plot at t={t:.3f} tstar={t*tscaling:.3f} (sec) c_int={c2[0]:.6f} a_int={alpha[0]:.6f} int={Lrefstar*hydride_interface[int(step/10)-1]} V={Lrefstar*(hydride_interface[int(step/10)-1])/(dt*tscaling)} c_end={c2[n2-1]} D_int={D[0]*D_factor:.4f} itn={iteration}"
    #     # )
    #     print(
    #                 f"Plot at t={t:.3f} tstar={t*tscaling:.3f} (sec) c_int={c2[0]:.6f} a_int={alpha[0]:.6f} c_end={c2[n2-1]} D_int={D[0]:.4f} itn={iteration},Gamma={Gamma}"
    #             )
    #     print(Gamma)
    #     # plt.plot(hydride_interface)
    #     # plt.pause(0.01)

    #     plt.plot(cs + reactK**(-1/2) * c2)
    #     plt.pause(0.001)

    #     np.savetxt(
    #         f"formal_work/data1D/inner_moving_c2{t:.3f}.dat",
    #         np.transpose(np.array([x2_nodes, c2])),
    #     )
    #     np.savetxt(
    #         f"formal_work/data1D/inner_moving_alpha{t:.3f}.dat",
    #         np.transpose(np.array([x2_nodes, alpha])),
    #     )
    #     # np.savetxt(
    #     #     f"1Dexamples/data/D_{t:.2f}.dat",
    #     #     np.transpose(np.array([x2_nodes, D])),
    #     # )
    # step += 1

