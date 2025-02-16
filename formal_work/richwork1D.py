
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
n3 = 51  # oxide
n2 = 101  # metal
# number of variables in the bulk
nv = 3  # [H] and [U] in the bulk plus diffusivity of mixture

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 1000 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 10 * 1.0e-7  # oxide domain lengthscale 10nm
D1star = 1.0e-13  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 1.49e-10  # cm^2/s, diffusion of H in U (room temp value)
D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 (room temp value)
N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 1.0e-5  # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castar = 1.0e-4  # mol/cm^3, surface value for [H] -- *ad-hoc* value

# fixed reference values to keep time/lengthscales consistent in comparisons
Lrefstar = 1.0e3 * 1.0e-7  # using 1um=1000nm here
Drefstar = D2star  # using the U value

# non-dimensional domains
L3 = L3star / Lrefstar  # oxide
L2 = L2star / Lrefstar  # bulk
# Choose computation coord Xmax such that x_from_X(Xmax)=L2
Xmax = sp.optimize.root(lambda X: x_from_X(X)[0] - L2, 5).x[0]
# exit if sanity check fails
assert abs(x_from_X(Xmax)[0] - L2) < 1.0e-8
# nodes in the oxide and bulk
x3_nodes = np.linspace(-L3, 0, n3)  # oxide, the interface is at x=0
X2_nodes = np.linspace(0, Xmax, n2)  # bulk nodes in comp. coordinate
# physical node locations
x2_nodes, x2d_nodes = x_from_X(X2_nodes)
# non-dimensional max = 10^4 sec
tmax = 3.0e4 * Drefstar / (Lrefstar * Lrefstar)

# non-dimensional variables
c3 = np.zeros(n3, dtype=np.float64)  # oxide
c2 = np.zeros(n2, dtype=np.float64)  # metal
alpha = np.ones(n2, dtype=np.float64)  # volume fraction of metal
D3 = D3star / Drefstar
D2 = D2star / Drefstar
D1 = D1star / Drefstar
cs = Csstar / Castar  # saturation value \in [0,1)
eps = Castar / N2star  # relative forcing measure
reactK = (
    1.0e6  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
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
# metric_file = open("1Dexamples/data/metric.dat", "w")

# step sizes
h3 = L3 / (n3 - 1)
h3s = h3 * h3
h2 = L2 / (n2 - 1)
h2s = h2 * h2
dt = 0.02
tol = 1.0e-8

# time step the system
t = 0
step = 1  # periodic output counter

while t < tmax:
    # evaluation at the new time step
    t += dt
    # store previous time step solution profiles
    c3_old = np.array(c3)
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

        # surface has [H]=c3=1 for node 0
        row += [0]
        col += [0]
        val += [1]
        b += [
            1 - np.exp(-t) - c3[0]
        ]  # ad-hoc smooth ramp up of forcing to 1 at surface

        # internal nodes: apply the diffusion equation in the oxide layer
        for j in range(1, n3 - 1):
            row += 3 * [j]
            col += [j - 1, j, j + 1]
            val += [D3 / h3s, -2 * D3 / h3s - 2 / dt, D3 / h3s]
            # residuals for a Crank-Nicolson method
            b += [
                2 * (c3[j] - c3_old[j]) / dt
                - D3 * (c3[j - 1] - 2 * c3[j] + c3[j + 1]) / h3s
                - D3 * (c3_old[j - 1] - 2 * c3_old[j] + c3_old[j + 1]) / h3s
            ]

        # interface c3 - c2 = 0 for node n3-1 in oxide : matched [H] at interface
        j = n3 - 1
        row += 2 * [j]
        col += [j, j + 1]
        val += [-1, 1]
        b += [c3[n3 - 1] - c2[0]]

        # internal nodes in metal/bulk region
        for j in range(0, n2):
            k = n3 + j * nv
            #
            ##################
            # DIFFUSION OF H #
            ##################
            #
            # at the first and last node we need BCns for the H diffusion
            if (j == 0) or (j == n2 - 1):
                if j == 0:
                    # at the first bulk node we impose the flux condition
                    row += 7 * [k]
                    col += [
                        k - 3,  # c3_{n3-3}
                        k - 2,  # c3_{n3-2}
                        k - 1,  # c3_{n3-1}
                        k,  # c2_0
                        k + 2,  #  D_0
                        k + nv,  # c2_1
                        k + 2 * nv,  # c2_2
                    ]
                    val += [
                        -D3 / (2 * h3),
                        +4 * D3 / (2 * h3),
                        -3 * D3 / (2 * h3),
                        -3 * D[0] / (2 * h2 * x2d_nodes[0]),
                        (-3 * c2[0] + 4 * c2[1] - c2[2]) / (2 * h2 * x2d_nodes[0]),
                        4 * D[0] / (2 * h2 * x2d_nodes[0]),
                        -D[0] / (2 * h2 * x2d_nodes[0]),
                    ]
                    # second order one sided derivatives
                    b += [
                        D3 * (3 * c3[n3 - 1] - 4 * c3[n3 - 2] + c3[n3 - 3]) / (2 * h3)
                        - D[0]
                        * (-3 * c2[0] + 4 * c2[1] - c2[2])
                        / (2 * h2 * x2d_nodes[0])  # deriv in terms of comp. coord.
                    ]
                if j == n2 - 1:
                    # at the last node we impose no H flux out of the domain
                    k = n3 + j * nv
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
                        - 3 * reactK * (3 * (c2[j] - cs) ** 2 * alpha[j]),  # c_j
                        -3 * reactK * (c2[j] - cs) ** 3,  # alpha_j
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
                        + 3 * reactK * ((c2[j] - cs) ** 3) * alpha[j]
                        # old values for Crank-Nicolson
                        - 0.5
                        * (D_old[j + 1] + D_old[j])
                        * (c2_old[j + 1] - c2_old[j])
                        / (h2s * x2d_nodes[j] * x2dph)
                        + 0.5
                        * (D_old[j] + D_old[j - 1])
                        * (c2_old[j] - c2_old[j - 1])
                        / (h2s * x2d_nodes[j] * x2dmh)
                        + 3 * reactK * ((c2_old[j] - cs) ** 3) * alpha_old[j]
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
                    -eps * reactK * (3 * (c2[j] - cs) ** 2 * alpha[j]),
                    -eps * reactK * (c2[j] - cs) ** 3 - 2 / dt,
                ]
                # residuals for a Crank-Nicolson method
                b += [
                    2 * (alpha[j] - alpha_old[j]) / dt
                    + eps * reactK * ((c2[j] - cs) ** 3) * alpha[j]
                    + eps * reactK * ((c2_old[j] - cs) ** 3) * alpha_old[j]
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
            val += [D2 - D1, -1]
            b += [D[j] - alpha[j] * D2 - (1 - alpha[j]) * D1]

        # solve the matrix problem for this iteration
       
        a = sp.sparse.coo_matrix((val, (row, col)), shape=(n3 + n2 * nv, n3 + n2 * nv))
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
        c3[0:n3] += x[0:n3]
        c2[0:n2] += x[n3 : n3 + n2 * nv : nv]
        alpha[0:n2] += x[n3 + 1 : n3 + n2 * nv : nv]
        D[0:n2] += x[n3 + 2 : n3 + n2 * nv : nv]

        if iteration > 10:
            print("Too many iterations")
            sys.exit()
    #
    # Solution for this time step has converged
    #

    # measure the interface value of [H] and position of [alpha]=0.01
    if step % 10 == 0:
        hydride_interface = 0
        # find positions where [UH3] concentration is 99% in the bulk
        interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - 0.01)
        roots = interpolated_alpha.roots(extrapolate=False)
        if len(roots) > 0:
            hydride_interface = roots[0]
        # OUTPUT:
        # non-dim time, dim time (sec), non-dim [H] @interface, non-dim [U] @interface, dim UH3 depth (nm)
        # metric_file.write(
        #     f"{t} {t*(Lrefstar**2)/(Drefstar)} {c2[0]} {alpha[0]} {hydride_interface*Lrefstar*1.e7}\n"
        # )
        # metric_file.flush()

    # save every 100 steps
    if step % 100 == 0:
        print(
            f"Plot at t={t} tstar={t*(Lrefstar**2)/(Drefstar)} (sec) c_int={c2[0]} c_end={c2[n2-1]} itn={iteration}"
        )
        np.savetxt(
            f"formal_work/data1D/k1e6/c3{t:.2f}.dat",
            np.transpose(np.array([x3_nodes, c3])),
        )
        np.savetxt(
            f"formal_work/data1D/k1e6/c2{t:.2f}.dat",
            np.transpose(np.array([x2_nodes, c2])),
        )
        np.savetxt(
            f"formal_work/data1D/k1e6/alpha{t:.2f}.dat",
            np.transpose(np.array([x2_nodes, alpha])),
        )
        # np.savetxt(
        #     f"1Dexamples/data/D_{t:.2f}.dat",
        #     np.transpose(np.array([x2_nodes, D])),
        # )
    step += 1

