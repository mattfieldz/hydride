
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
n2 = 101 # metal
m = 40 # width nodes

# number of variables in the bulk
nv = 3  # [H] and [U] in the bulk plus diffusivity of mixture

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 1000 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 10 * 1.0e-7  # oxide domain lengthscale 10nm
Lmstar = 1000 * 1.0e-7 # 100um
D1star = 1.0e-15  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
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
Lm = Lmstar / L2star # width

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

np.savetxt(
            f"formal_work/data2D/k1e6/single/x2_nodes.dat",
            x2_nodes,
        )
np.savetxt(
            f"formal_work/data2D/k1e6/single/x3_nodes.dat",
            x3_nodes,
        )
# non-dimensional variables
c3 = np.zeros((n3,m), dtype=np.float64)  # oxide
c2 = np.zeros((n2,m), dtype=np.float64)  # metal
alpha = np.ones((n2,m), dtype=np.float64)  # volume fraction of metal
D3 = D3star / Drefstar
D2 = D2star / Drefstar
D1 = D1star / Drefstar
cs = Csstar / Castar  # saturation value \in [0,1)
eps = Castar / N2star  # relative forcing measure
reactK = (
    1e6  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
)
D = D2 * np.ones((n2,m), dtype=np.float64)  # mixture diffusion in the bulk, initially =D2

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
hl = Lm / (m-1)
hls = hl * hl

dt = 0.02
tol = 1.0e-8

# time step the system
t = 0
step = 1  # periodic output counter

while t < tmax:
    # evaluation at the new time step
    t += dt
    # store previous time step solution profiles
    c3_old = np.copy(c3)
    c2_old = np.copy(c2)
    alpha_old = np.copy(alpha)
    D_old = np.copy(D)

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
        count = 0
        for m_i in range(0,m):
            m_s = (n3 + nv * n2)*m_i # current node
            m_j = (n3 + nv * n2) # spacing
            row += [0+m_s]
            col += [0+m_s]
            val += [1]
            b += [
                1 - np.exp(-t) - c3[0,m_i]
            ]  # ad-hoc smooth ramp up of forcing to 1 at surface
            
            # internal nodes: apply the diffusion equation in the oxide layer
            
            
            
                        
                    
                
                        
                    
                        
                        
            # surface has [H]=c3=1 for node 0
        
    
            for j in range(1, n3 - 1):
                if m_i == m-1:
                    row += 3 * [j + m_s]
                    col += [(j-2*m_j) + m_s, (j-m_j) + m_s, (j)+m_s]
                    val += [-1, 4, -3]
                    b += [c3[j,m_i-2] + 3*c3[j,m_i] - 4*c3[j,m_i-1]]
                elif m_i == 0:
                    row += 3 * [j + m_s]
                    col += [(j+2*m_j) + m_s, (j+m_j) + m_s, (j)+m_s]
                    val += [-1, 4, -3]
                    b += [c3[j,m_i+2] + 3*c3[j,m_i] - 4*c3[j,m_i+1]]
                elif (j == 1) and (18 < m_i) and (m_i < 21):
                    
                    row += [j + m_s]
                    col += [j + m_s]
                    val += [1]
                    b += [1 - np.exp(-t) - c3[j,m_i]]
                
                else:
                    
                    row += 5 * [j+m_s]
                    col += [j - 1 + m_s, j + m_s, j + 1 + m_s,j + m_s + m_j , j + m_s - m_j]
                    val += [D3 / h3s, -2 * D3 / h3s - 2 * D3 / hls - 2 / dt, D3 / h3s, D3 / hls, D3/hls]
                    # residuals for a Crank-Nicolson method
                    b += [
                        2 * (c3[j,m_i] - c3_old[j,m_i]) / dt
                        - D3 * (c3[j - 1,m_i] - 2 * c3[j,m_i] + c3[j + 1,m_i]) / h3s
                        - D3 * (c3_old[j - 1,m_i] - 2 * c3_old[j,m_i] + c3_old[j + 1,m_i]) / h3s
                        - D3 * (c3[j,m_i-1] - 2 * c3[j,m_i] + c3[j,m_i+1]) / hls
                        - D3 * (c3_old[j,m_i-1] - 2 * c3_old[j,m_i] + c3_old[j,m_i+1]) / hls
                    ]
                
            # interface c3 - c2 = 0 for node n3-1 in oxide : matched [H] at interface
            j = n3 - 1
            row += 2 * [j + m_s]
            col += [j + m_s, j + 1 + m_s]
            val += [-1, 1]
            b += [c3[n3 - 1,m_i] - c2[0,m_i]]
            # internal nodes in metal/bulk region
            for j in range(0, n2):
                k = n3 + j * nv + m_s
                
                
                #
                ##################
                # DIFFUSION OF H #
                ##################
                #

                if m_i == 0:
                    
                    row += 3 * [k]
                    col += [(k+2*m_j), (k+m_j), (k)]
                    val += [-1, 4, -3]
                    b += [c2[j,m_i+2] + 3*c2[j,m_i] - 4*c2[j,m_i+1]]
                elif m_i == m-1:
                    
                
                    row += 3 * [k]
                    col += [(k-2*m_j), (k-m_j), (k)]
                    val += [-1, 4, -3]
                    b += [c2[j,m_i-2] + 3*c2[j,m_i] - 4*c2[j,m_i-1]]
                    
                else:    
                
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
                                -3 * D[0,m_i] / (2 * h2 * x2d_nodes[0]),
                                (-3 * c2[0,m_i] + 4 * c2[1,m_i] - c2[2,m_i]) / (2 * h2 * x2d_nodes[0]),
                                4 * D[0,m_i] / (2 * h2 * x2d_nodes[0]),
                                -D[0,m_i] / (2 * h2 * x2d_nodes[0]),
                            ]
                            # second order one sided derivatives
                            b += [
                                D3 * (3 * c3[n3 - 1,m_i] - 4 * c3[n3 - 2,m_i] + c3[n3 - 3,m_i]) / (2 * h3)
                                - D[0,m_i]
                                * (-3 * c2[0,m_i] + 4 * c2[1,m_i] - c2[2,m_i])
                                / (2 * h2 * x2d_nodes[0])  # deriv in terms of comp. coord.
                            ]
                            
                            
                        if j == n2 - 1:
                            # at the last node we impose no H flux out of the domain
                            
                            row += 3 * [k]
                            col += [k - 2 * nv, k - nv, k]
                            val += [-1, 4, -3]
                            # second order differencing
                            b += [(3 * c2[j,m_i] - 4 * c2[j - 1,m_i] + c2[j - 2,m_i])]
                            
                    else:  # we impose the diffusion equation at internal nodes
                        # these values appear in the diffusion w.r.t. comp. coord.
                        x2dph = 0.5 * (x2d_nodes[j + 1] + x2d_nodes[j])
                        x2dmh = 0.5 * (x2d_nodes[j] + x2d_nodes[j - 1])
                        if c2[j,m_i] > cs:
                            # [H] diffusion equation, reaction is active
                            # with local duffisivity D[j] = D2*alpha[j] + D1*(1-alpha[j])
                            row += 11 * [k]
                            
                            col += [
                                k - nv,  # c_{j-1}
                                k - nv + 2,  # D_{j-1}
                                k,  # c_j
                                k + 1,  # alpha_j
                                k + 2,  # D_j
                                k + nv,  # c_{j+1}
                                k + nv + 2,  # D_{j+1}
                                k + m_j, # right node
                                k - m_j, # left node
                                k + m_j + 2,
                                k - m_j + 2
                            ]
                            val += [
                                0.5
                                * (D[j,m_i] + D[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh),  # c_{j-1}
                                -0.5
                                * (c2[j,m_i] - c2[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh),  # D_{j-1}
                                -0.5
                                * ((D[j + 1,m_i] + D[j,m_i]) / x2dph + (D[j,m_i] + D[j - 1,m_i]) / x2dmh)
                                / (h2s * x2d_nodes[j])
                                -0.5
                                * ((D[j,m_i+1] + D[j,m_i]) + (D[j,m_i] + D[j,m_i-1]))
                                / (hls) # left ride contributions
                                - 2 / dt
                                - 3 * reactK * (3 * (c2[j,m_i] - cs) ** 2 * alpha[j,m_i]),  # c_j
                                -3 * reactK * (c2[j,m_i] - cs) ** 3,  # alpha_j
                                +0.5 * (c2[j + 1,m_i] - c2[j,m_i]) / (h2s * x2d_nodes[j] * x2dph)
                                - 0.5
                                * (c2[j,m_i] - c2[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh) 
                                +0.5 * (c2[j,m_i+1] - c2[j,m_i]) / (hls) # left ride contributions
                                - 0.5
                                * (c2[j,m_i] - c2[j,m_i-1])
                                / (hls),  # D_j
                                0.5
                                * (D[j + 1,m_i] + D[j,m_i])
                                / (h2s * x2d_nodes[j] * x2dph),  # c_{j+1}
                                +0.5
                                * (c2[j + 1,m_i] - c2[j,m_i])
                                / (h2s * x2d_nodes[j] * x2dph),  # D_{j+1}
                                0.5
                                * (D[j,m_i+1] + D[j,m_i])
                                / (hls), # right node
                                0.5
                                * (D[j,m_i] + D[j,m_i-1])
                                / (hls), # left node
                                +0.5
                                * (c2[j, m_i+1] - c2[j,m_i])
                                / (hls),  # D_{j+1} right node
                                +0.5
                                * (c2[j,m_i] - c2[j,m_i-1])
                                / (hls)  # D_{j+1} left node
                                
                            ]
                            # residuals for Crank-Nicolson method
                            b += [
                                2 * (c2[j,m_i] - c2_old[j,m_i]) / dt
                                - 0.5
                                * (D[j + 1,m_i] + D[j,m_i])
                                * (c2[j + 1,m_i] - c2[j,m_i])
                                / (h2s * x2d_nodes[j] * x2dph)
                                + 0.5
                                * (D[j,m_i] + D[j - 1,m_i])
                                * (c2[j,m_i] - c2[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh)

                                - 0.5
                                * (D[j,m_i+1] + D[j,m_i])
                                * (c2[j,m_i+1] - c2[j,m_i])
                                / (hls)
                                + 0.5
                                * (D[j,m_i] + D[j,m_i-1])
                                * (c2[j,m_i] - c2[j,m_i-1])
                                / (hls)

                                + 3 * reactK * ((c2[j,m_i] - cs) ** 3) * alpha[j,m_i]
                                # old values for Crank-Nicolson
                                - 0.5
                                * (D_old[j + 1,m_i] + D_old[j,m_i])
                                * (c2_old[j + 1,m_i] - c2_old[j,m_i])
                                / (h2s * x2d_nodes[j] * x2dph)
                                + 0.5
                                * (D_old[j,m_i] + D_old[j - 1,m_i])
                                * (c2_old[j,m_i] - c2_old[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh)

                                - 0.5
                                * (D_old[j,m_i+1] + D_old[j,m_i])
                                * (c2_old[j,m_i+1] - c2_old[j,m_i])
                                / (hls)
                                + 0.5
                                * (D_old[j,m_i] + D_old[j,m_i-1])
                                * (c2_old[j,m_i] - c2_old[j,m_i-1])
                                / (hls)

                                + 3 * reactK * ((c2_old[j,m_i] - cs) ** 3) * alpha_old[j,m_i]
                            ]
                            
                        else:
                            # [H] equation, no reaction
                            row += 10 * [k]
                            
                            col += [
                                k - nv,  # c_{j-1}
                                k - nv + 2,  # D_{j-1}
                                k,  # c_j
                                k + 2,  # D_j
                                k + nv,  # c_{j+1}
                                k + nv + 2,  # D_{j+1}
                                k + m_j, # right node
                                k - m_j, # left node
                                k + m_j + 2,
                                k - m_j + 2
                            ]
                            val += [
                                0.5
                                * (D[j,m_i] + D[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh),  # c_{j-1}
                                -0.5
                                * (c2[j,m_i] - c2[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh),  # D_{j-1}
                                -0.5
                                * ((D[j + 1,m_i] + D[j,m_i]) / x2dph + (D[j,m_i] + D[j - 1,m_i]) / x2dmh)
                                / (h2s * x2d_nodes[j])
                                -0.5
                                * ((D[j,m_i+1] + D[j,m_i]) + (D[j,m_i] + D[j,m_i-1]))
                                / (hls)
                                - 2 / dt,  # c_j
                                +0.5 * (c2[j + 1,m_i] - c2[j,m_i]) / (h2s * x2d_nodes[j] * x2dph)
                                - 0.5
                                * (c2[j,m_i] - c2[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh)
                                +0.5 * (c2[j,m_i+1] - c2[j,m_i]) / (hls)
                                - 0.5
                                * (c2[j,m_i] - c2[j,m_i-1])
                                / (hls),  # D_j
                                0.5
                                * (D[j + 1,m_i] + D[j,m_i])
                                / (h2s * x2d_nodes[j] * x2dph),  # c_{j+1}
                                +0.5
                                * (c2[j + 1,m_i] - c2[j,m_i])
                                / (h2s * x2d_nodes[j] * x2dph),  # D_{j+1}

                                0.5
                                * (D[j,m_i+1] + D[j,m_i])
                                / (hls),  # c_{j+1}

                                0.5
                                * (D[j,m_i] + D[j,m_i-1])
                                / (hls),  # c_{j-1}

                                +0.5
                                * (c2[j,m_i+1] - c2[j,m_i])
                                / (hls),  # D_{j+1}
                                
                                -0.5
                                * (c2[j,m_i] - c2[j,m_i-1])
                                / (hls)

                            ]
                            # residuals for Crank-Nicolson method
                            b += [
                                2 * (c2[j,m_i] - c2_old[j,m_i]) / dt
                                - 0.5
                                * (D[j + 1,m_i] + D[j,m_i])
                                * (c2[j + 1,m_i] - c2[j,m_i])
                                / (h2s * x2d_nodes[j] * x2dph)
                                + 0.5
                                * (D[j,m_i] + D[j - 1,m_i])
                                * (c2[j,m_i] - c2[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh)
                                - 0.5
                                * (D[j,m_i+1] + D[j,m_i])
                                * (c2[j,m_i+1] - c2[j,m_i])
                                / (hls)
                                + 0.5
                                * (D[j,m_i] + D[j,m_i-1])
                                * (c2[j,m_i] - c2[j,m_i-1])
                                / (hls)
                                # old time values for Crank-Nicolson
                                - 0.5
                                * (D_old[j + 1,m_i] + D_old[j,m_i])
                                * (c2_old[j + 1,m_i] - c2_old[j,m_i])
                                / (h2s * x2d_nodes[j] * x2dph)
                                + 0.5
                                * (D_old[j,m_i] + D_old[j - 1,m_i])
                                * (c2_old[j,m_i] - c2_old[j - 1,m_i])
                                / (h2s * x2d_nodes[j] * x2dmh)
                                - 0.5
                                * (D_old[j,m_i+1] + D_old[j,m_i])
                                * (c2_old[j,m_i+1] - c2_old[j,m_i])
                                / (hls)
                                + 0.5
                                * (D_old[j,m_i] + D_old[j,m_i-1])
                                * (c2_old[j,m_i] - c2_old[j,m_i-1])
                                / (hls)
                                
                            ]
                        
            #              
            ########################
            # CONSUMPTION OF U EQN #
            ########################
            #
        
                if c2[j,m_i] > cs:
                    # [alpha] equation, reaction is active
                    #
                    #
                    row += 2 * [k + 1]
                    col += [k, k + 1]  # c_{j},alpha_j
                    val += [
                        -eps * reactK * (3 * (c2[j,m_i] - cs) ** 2 * alpha[j,m_i]),
                        -eps * reactK * (c2[j,m_i] - cs) ** 3 - 2 / dt,
                    ]
                    # residuals for a Crank-Nicolson method
                    b += [
                        2 * (alpha[j,m_i] - alpha_old[j,m_i]) / dt
                        + eps * reactK * ((c2[j,m_i] - cs) ** 3) * alpha[j,m_i]
                        + eps * reactK * ((c2_old[j,m_i] - cs) ** 3) * alpha_old[j,m_i]
                    ]
                    
                else:
                    # [alpha] equation, no reaction terms
                    #
                    #
                    row += [k + 1]
                    col += [k + 1]  # alpha_j
                    
                    val += [-1 / dt]
                    b += [(alpha[j,m_i] - alpha_old[j,m_i]) / dt]
                    
                #
                ##############################
                # SET LOCAL DIFFUSION COEFFT #
                ##############################
                #
                row += 2 * [k + 2]
                col += [k + 1, k + 2]  # alpha_j,D_j
                
                val += [D2 - D1, -1]
                b += [D[j,m_i] - alpha[j,m_i] * D2 - (1 - alpha[j,m_i]) * D1]
                

                    
        # solve the matrix problem for this iteration
        
       
        a = sp.sparse.coo_matrix((val, (row, col)), shape=((n3 + n2 * nv) * m, (n3 + n2 * nv) * m))
        # system = petsc.PETScSparseLinearSystem(a, b)
        # x = system.linear_solve()
        # plt.imshow(a.toarray().astype(bool))
        # plt.show()
        # print('det',np.linalg.det(a.toarray()))
        x = sp.sparse.linalg.spsolve(a.tocsr(),b)
        # residual is the maximum correction value
        residual = sp.linalg.norm(x, ord=np.inf)

        
        
        # print(
        #    f"iteration = {iteration} residual = {residual} ||b||_inf = {sp.linalg.norm(b, ord=np.inf)}"
        # )
        # add the corrections to the current guess
        x_reshaped = x.reshape((n3+nv*n2,m),order='F')
        c3[0:n3] += x_reshaped[0:n3,:]
        c2[0:n2] += x_reshaped[n3 : n3 + n2 * nv : nv,:]
        alpha[0:n2] += x_reshaped[n3 + 1 : n3 + n2 * nv : nv,:]
        D[0:n2] += x_reshaped[n3 + 2 : n3 + n2 * nv : nv,:]

        if iteration > 10:
            print("Too many iterations")
            sys.exit()
        
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
        
        print('iterations : ',iteration)

        np.savetxt(
            f"formal_work/data2D/k1e6/single/c3{t:.2f}.dat",
            c3,
        )
        np.savetxt(
            f"formal_work/data2D/k1e6/single/c2{t:.2f}.dat",
            c2,
        )
        np.savetxt(
            f"formal_work/data2D/k1e6/single/alpha{t:.2f}.dat",
            alpha,
        )
        # np.savetxt(
        #     f"1Dexamples/data/D_{t:.2f}.dat",
        #     np.transpose(np.array([x2_nodes, D])),
        # )
    step += 1

