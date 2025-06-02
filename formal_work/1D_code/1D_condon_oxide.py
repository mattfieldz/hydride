
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

n2 = 4001  # metal
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

D_factor = 1
D_large = 1
Drefstar = D_factor * D2star # using the U value

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

c2 = np.zeros(n2, dtype=np.float64)  # metal
alpha = np.ones(n2, dtype=np.float64)  # volume fraction of metal
D3 = D3star / Drefstar
D2 = D2star / Drefstar
# D1 = D1star / Drefstar
D1 = D2
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

h2 = L2 / (Xmax*(n2 - 1))
h2s = h2 * h2
dt_standard = 0.00000003
dt = dt_standard
tol = 1.0e-8
threshold_alpha = 0.98



print(x2_nodes[200])

# time step the system
t = 0
step = 1  # periodic output counter

hydride_interface = []

t_0 = np.zeros(n2)

while t < tmax:
    # evaluation at the new time step
    t += dt
    # store previous time step solution profiles
    
    c2_old = np.array(c2)
    alpha_old = np.array(alpha)
    D_old = np.array(D)

    residual = 1
    iteration = 0
    # Newton iteration loop
    while residual > tol:
        iteration += 1
        # construct the sparse matrix
        row = []
        col = []
        val = []
        b = []

        # internal nodes in metal/bulk region
        for j in range(0, n2):
            k = j * nv
            #
            ##################
            # DIFFUSION OF H #
            ##################
            #
            # at the first and last node we need BCns for the H diffusion
            if (j == 0) or (j == n2 - 1):
                if j == 0:
                    # at the first bulk node we impose the flux condition
                    row += 3 * [k]
                    col += [
                        
                        k,  # c2_0
                         #  D_0
                        k + nv,  # c2_1
                        k + 2 * nv,  # c2_2
                    ]
                    val += [
                        -3 * D[0] / (2 * h2 * x2d_nodes[0]),
                        4 * D[0] / (2 * h2 * x2d_nodes[0]),
                        -D[0] / (2 * h2 * x2d_nodes[0]),
                    ]
                    # second order one sided derivatives
                    b += [
                        D3 * (c2[j]-1) / L3
                        - D[0]
                        * (-3 * c2[0] + 4 * c2[1] - c2[2])
                        / (2 * h2 * x2d_nodes[0])  # deriv in terms of comp. coord.
                    ]
                if j == n2 - 1:
                    # at the last node we impose no H flux out of the domain
                    k = j * nv
                    row += 3 * [k]
                    col += [k - 2 * nv, k - nv, k]
                    val += [-1, 4, -3]
                    # second order differencing
                    b += [(3 * c2[j] - 4 * c2[j - 1] + c2[j - 2])]
            else:  # we impose the diffusion equation at internal nodes
                # these values appear in the diffusion w.r.t. comp. coord.



                x2dph = 0.5 * (x2d_nodes[j + 1] + x2d_nodes[j])
                x2dmh = 0.5 * (x2d_nodes[j] + x2d_nodes[j - 1])

                
                
                
                if c2[j] > cs:
                    # [H] diffusion equation, reaction is active
                    # with local duffisivity D[j] = D2*alpha[j] + D1*(1-alpha[j])
                    row += 7 * [k]
                    col += [
                        k - nv,  # c_{j-1}
                        k - nv + 2,  # D_{j-1}
                        k,  # c_j
                        k + 1,  # alpha_j
                        k + 2,  # D_j
                        k + nv,  # c_{j+1}
                        k + nv + 2,  # D_{j+1}
                    ]
                    val += [
                        0.5
                        * (D[j] + D[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh),  # c_{j-1}
                        -0.5
                        * (c2[j] - c2[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh),  # D_{j-1}
                        -0.5
                        * ((D[j + 1] + D[j]) / x2dph + (D[j] + D[j - 1]) / x2dmh)
                        / (h2s * x2d_nodes[j])
                        - 2 / dt
                        - 3 * reactK * (m * (c2[j] - cs) ** (m-1) * alpha[j]),  # c_j
                        -3 * reactK * (c2[j] - cs) ** m,  # alpha_j
                        +0.5 * (c2[j + 1] - c2[j]) / (h2s * x2d_nodes[j] * x2dph)
                        - 0.5
                        * (c2[j] - c2[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh),  # D_j
                        0.5
                        * (D[j + 1] + D[j])
                        / (h2s * x2d_nodes[j] * x2dph),  # c_{j+1}
                        +0.5
                        * (c2[j + 1] - c2[j])
                        / (h2s * x2d_nodes[j] * x2dph),  # D_{j+1}
                    ]
                    # residuals for Crank-Nicolson method
                    b += [
                        2 * (c2[j] - c2_old[j]) / dt
                        - 0.5
                        * (D[j + 1] + D[j])
                        * (c2[j + 1] - c2[j])
                        / (h2s * x2d_nodes[j] * x2dph)
                        + 0.5
                        * (D[j] + D[j - 1])
                        * (c2[j] - c2[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh)
                        + 3 * reactK * ((c2[j] - cs) ** m) * alpha[j]
                        # old values for Crank-Nicolson
                        - 0.5
                        * (D_old[j + 1] + D_old[j])
                        * (c2_old[j + 1] - c2_old[j])
                        / (h2s * x2d_nodes[j] * x2dph)
                        + 0.5
                        * (D_old[j] + D_old[j - 1])
                        * (c2_old[j] - c2_old[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh)
                        + 3 * reactK * ((c2_old[j] - cs) ** m) * alpha_old[j]
                    ]
                else:
                    # [H] equation, no reaction
                    row += 6 * [k]
                    col += [
                        k - nv,  # c_{j-1}
                        k - nv + 2,  # D_{j-1}
                        k,  # c_j
                        k + 2,  # D_j
                        k + nv,  # c_{j+1}
                        k + nv + 2,  # D_{j+1}
                    ]
                    val += [
                        0.5
                        * (D[j] + D[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh),  # c_{j-1}
                        -0.5
                        * (c2[j] - c2[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh),  # D_{j-1}
                        -0.5
                        * ((D[j + 1] + D[j]) / x2dph + (D[j] + D[j - 1]) / x2dmh)
                        / (h2s * x2d_nodes[j])
                        - 2 / dt,  # c_j
                        +0.5 * (c2[j + 1] - c2[j]) / (h2s * x2d_nodes[j] * x2dph)
                        - 0.5
                        * (c2[j] - c2[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh),  # D_j
                        0.5
                        * (D[j + 1] + D[j])
                        / (h2s * x2d_nodes[j] * x2dph),  # c_{j+1}
                        +0.5
                        * (c2[j + 1] - c2[j])
                        / (h2s * x2d_nodes[j] * x2dph),  # D_{j+1}
                    ]
                    # residuals for Crank-Nicolson method
                    b += [
                        2 * (c2[j] - c2_old[j]) / dt
                        - 0.5
                        * (D[j + 1] + D[j])
                        * (c2[j + 1] - c2[j])
                        / (h2s * x2d_nodes[j] * x2dph)
                        + 0.5
                        * (D[j] + D[j - 1])
                        * (c2[j] - c2[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh)
                        # old time values for Crank-Nicolson
                        - 0.5
                        * (D_old[j + 1] + D_old[j])
                        * (c2_old[j + 1] - c2_old[j])
                        / (h2s * x2d_nodes[j] * x2dph)
                        + 0.5
                        * (D_old[j] + D_old[j - 1])
                        * (c2_old[j] - c2_old[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh)
                    ]
            #
            ########################
            # CONSUMPTION OF U EQN #
            ########################
            #
            if c2[j] > cs:
                # [alpha] equation, reaction is active
                #
                #
                row += 2 * [k + 1]
                col += [k, k + 1]  # c_{j},alpha_j
                val += [
                    -eps * reactK * (m * (c2[j] - cs) ** (m-1) * alpha[j]),
                    -eps * reactK * (c2[j] - cs) ** m - 2 / dt,
                ]
                # residuals for a Crank-Nicolson method
                b += [
                    2 * (alpha[j] - alpha_old[j]) / dt
                    + eps * reactK * ((c2[j] - cs) ** m) * alpha[j]
                    + eps * reactK * ((c2_old[j] - cs) ** m) * alpha_old[j]
                ]
            else:
                # [alpha] equation, no reaction terms
                #
                #
                row += [k + 1]
                col += [k + 1]  # alpha_j
                val += [-1 / dt]
                b += [(alpha[j] - alpha_old[j]) / dt]
            #
            ##############################
            # SET LOCAL DIFFUSION COEFFT #
            ##############################
            #
            row += 2 * [k + 2]
            col += [k + 1, k + 2]  # alpha_j,D_j
            val += [D2 - D_large, -1]
            b += [D[j] - alpha[j] * D2 - (1 - alpha[j]) * D_large]
            
            # if alpha[j] > threshold_alpha:
            #     row += [k + 2]
            #     col += [k + 2]  # alpha_j,D_j
            #     val += [-1]
            #     b += [D[j] - D2]
            # else:
            #     if t_0[j] == 0:
            #         t_0[j] = t
                
            #     row += [k + 2]
            #     col += [k + 2]  # alpha_j,D_j
            #     val += [-1]
            #     # b += [D[j] - D2 * (1 + ((t-t_0[j])/15)**3)]
            #     if alpha[j] > 0.8:
            #         b += [D[j] - (D2 + D_large / (0.8-threshold_alpha)**3 * (alpha[j]-threshold_alpha)**3)]
            #     else:
            #         b += [D[j] - D_large]
        # solve the matrix problem for this iteration
       
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
        D[0:n2] += x[2 : n2 * nv : nv]

        # if iteration > 5:
        #     # print("Too many iterations")
        #     dt = dt / 2
        # elif iteration < 3:
        #     # print("Too few iterations")
        #     dt = dt * 2
        # else:
        #     dt = dt_standard


        if iteration > 10:
            print('too many')
            print(residual)
    #   
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
    if step % 100 == 0:
        # print(
        #     f"Plot at t={t:.3f} tstar={t*tscaling:.3f} (sec) c_int={c2[0]:.6f} a_int={alpha[0]:.6f} int={Lrefstar*hydride_interface[int(step/10)-1]} V={Lrefstar*(hydride_interface[int(step/10)-1])/(dt*tscaling)} c_end={c2[n2-1]} D_int={D[0]*D_factor:.4f} itn={iteration}"
        # )
        print(
                    f"Plot at t={t:.7f} tstar={t*tscaling:.3f} (sec) c_int={c2[0]:.6f} a_int={alpha[0]:.6f} c_end={c2[n2-1]} D_int={D[0]*D_factor:.4f} itn={iteration}"
                )
        # plt.plot(hydride_interface)
        # plt.pause(0.01)

        np.savetxt(
            f"formal_work/data1D/c2_condon_oxide_no_threshold_D1e0{t:.7f}.dat",
            np.transpose(np.array([x2_nodes, c2])),
        )
        np.savetxt(
            f"formal_work/data1D/alpha_condon_oxide_no_threshold_D1e0{t:.7f}.dat",
            np.transpose(np.array([x2_nodes, alpha])),
        )
        # np.savetxt(
        #     f"1Dexamples/data/D_{t:.2f}.dat",
        #     np.transpose(np.array([x2_nodes, D])),
        # )
    step += 1

