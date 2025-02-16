import numpy as np
import scipy as sp
import sys 
# Crank-Nicholson approach to solving reaction diffusion equation in 1 dimension iterating

# Function as in Rich's work for rescaling the coordinates in the metal region to resolve at the interface
def x_from_X(X):
    """Defines the (non-uniformly distributed) physical coordinate 'x'
    in terms of the (uniformly distributed) computational coordinate
    'X'. Returns: not only x(X) but also x'(X) and x''(X).

    """
    BX = 0.01  # smaller BX leads to more NON-uniformity
    # +0*X below prevents numpy from returning a scalar when acting
    # on an iterable container for X
    return (X + BX) * (X + BX) - BX * BX, 2 * (X + BX)

# oxide-mesh
n3 = 51
# metal-mesh
n2 = 1001
# Value to skip between rows of matrix
ej = 2





# Dimensional quantities

# Diffusion coefficients in cm^2/s
D1star = 1e-13
D2star = 1e-10
D3star = 1e-12

# Length scale for oxide and metal in cm

L3star = 10 * 1e-7 # 10nm
L2star = 1000 * 1e-7 # 1000nm

# Reference concentration values molcm^-3

Castar = 1e-4 # surface hydrogen concentration
Csstar = 1e-5 # solubility limit ~ 0.1

# Number density of uranium metal molcm^-3

N2star = 8.01e-2 

# Dimensional reaction coefficient
kstar = 1e15


# Reference length scale and diffusivity

Dref = D2star
Lref = L2star

D1 = D1star/D2star
D2 = 1
D3 = D3star/D2star

L3 = L3star/L2star
L2 = 1

# array of hydrogen concentrations
c3 = np.zeros(n3)
c2 = np.zeros(n2)

# array of uranium concentrations
alpha = np.ones(n2)

# array of diffusivities
D = D2 * np.ones(n2)

# Choose computation coord Xmax such that x_from_X(Xmax)=L2
Xmax = sp.optimize.root(lambda X: x_from_X(X)[0] - L2, 5).x[0]
# exit if sanity check fails
assert abs(x_from_X(Xmax)[0] - L2) < 1.0e-8
# nodes in the oxide and bulk
x3_nodes = np.linspace(-L3, 0, n3)  # oxide, the interface is at x=0
X2_nodes = np.linspace(0, Xmax, n2)  # bulk nodes in comp. coordinate
# physical node locations
x2_nodes, x2d_nodes = x_from_X(X2_nodes)




cs = Csstar / Castar # solubility limit
eps = Castar / N2star

k = 1e5


# maximum time
t_max = 20

# step size
h3 = L3 / (n3 - 1)
h3s = h3 * h3
h2 = L2 / (n2 - 1)
h2s = h2 * h2
dt = 0.01
tolerance = 1.0e-8

# time 
t = 0
# counter 
step = 1  

while t < t_max:
    # old values for crank nicholson
    t += dt
    c2_o = np.copy(c2)
    c3_o = np.copy(c3)
    alpha_o = np.copy(alpha)
    D_o = np.copy(D)
    res = 1
    iterations = 0

    while res > tolerance:
        
        iterations += 1
        # Computes diffusivities
        
        row = []
        col = []
        val = []
        b = []

        # gas-oxide interface
        row += [0]
        col += [0]
        val += [1]
        b += [1 - np.exp(-t) - c3[0]] # c + c_tilde = 1 at surface
        
        # main bulk of the oxide no reaction

        for i in range(1,n3-1):
            row += 3 * [i]
            col += [i-1,i,i+1]
            val += [D3/h3s,-2 * D3/h3s - 2/dt, D3/h3s] # finite difference quotient
            b += [2/dt * (c3[i]-c3_o[i]) - 
                D3/h3s * (c3[i-1] - 2 * c3[i] + c3)[i+1] -
                D3/h3s * (c3_o[i-1] - 2 * c3_o[i] + c3_o)[i+1]] # crank nicholson residuals
        
        # interface matching concentrations
        i = n3-1
        row += 2 * [i]
        col += [i,i+1]
        val += [-1,1]
        b += [c3[n3-1] - c2[0]]

        # interface matching gradients (Fick's second law)
        i = n3 

        
        row += 6 * [i]
        col += [i-3,i-2,i-1,i,i+ej,i+2*ej] # Either side of the boundary stencil
        val += [-D3 / (2 * h3),
                +4 * D3 / (2 * h3), # oxide boundary 
                -3 * D3 / (2 * h3),
                
                -3 * D[0] / (2 * h2 * x2d_nodes[0]),        
                4 * D[0] / (2 * h2 * x2d_nodes[0]), # metal/hydride boundary
                -D[0] / (2 * h2 * x2d_nodes[0])]
        b += [ D3 * (3 * c3[n3 - 1] - 4 * c3[n3 - 2] + c3[n3 - 3]) / (2 * h3)
                - D[0] * (-3 * c2[0] + 4 * c2[1] - c2[2]) / (2 * h2 * x2d_nodes[0])]
        
        for i in range(0,n2):
            j = n3 + ej * i # scales up the number of new matrix columns/rows noting that we need extra entries for alpha terms
            
            if (i != 0) and (i != n2-1):
                if c2[i] < cs:  # no reaction
                    x2dph = 0.5 * (x2d_nodes[i + 1] + x2d_nodes[i])
                    x2dmh = 0.5 * (x2d_nodes[i] + x2d_nodes[i - 1])
                    
                    row += 3 * [j]
                    col += [
                        j - ej,  # c_{j-1}
                        j,  # c_j
                        j + ej,  # c_{j+1}
                    ]
                    val += [
                        0.5
                        * (D[i] + D[i - 1])
                        / (h2s * x2d_nodes[i] * x2dmh),  # c_{j-1}

                        -0.5
                        * ((D[i + 1] + D[i]) / x2dph + (D[i] + D[i - 1]) / x2dmh)
                        / (h2s * x2d_nodes[i])
                        - 2 / dt,  # c_j

                        0.5
                        * (D[i + 1] + D[i])
                        / (h2s * x2d_nodes[i] * x2dph),  # c_{j+1}

                        ]
                        # residuals for Crank-Nicolson method
                    b += [
                        2 * (c2[i] - c2[i]) / dt
                        - 0.5
                        * (D[i + 1] + D[i])
                        * (c2[i + 1] - c2[i])
                        / (h2s * x2d_nodes[i] * x2dph)
                        + 0.5
                        * (D[i] + D[i - 1])
                        * (c2[i] - c2[i - 1])
                        / (h2s * x2d_nodes[i] * x2dmh)
                        # old time values for Crank-Nicolson
                        - 0.5
                        * (D_o[i + 1] + D_o[i])
                        * (c2_o[i + 1] - c2_o[i])
                        / (h2s * x2d_nodes[i] * x2dph)
                        + 0.5
                        * (D_o[i] + D_o[i - 1])
                        * (c2_o[i] - c2_o[i - 1])
                        / (h2s * x2d_nodes[i] * x2dmh)
                    ]
                else:
                    # [H] diffusion equation, reaction is active
                        # with local duffisivity D[i] = D2*alpha[i] + D1*(1-alpha[i])
                    row += 4 * [j]
                    col += [
                        j - ej,  # c_{j-1}
                        
                        j,  # c_j
                        j + 1,  # alpha_j
                        
                        j + ej,  # c_{j+1}
                        
                    ]
                    val += [
                        0.5
                        * (D[i] + D[i - 1])
                        / (h2s * x2d_nodes[i] * x2dmh),  # c_{j-1}
            
                        -0.5
                        * ((D[i + 1] + D[i]) / x2dph + (D[i] + D[i - 1]) / x2dmh)
                        / (h2s * x2d_nodes[i])
                        - 2 / dt
                        - 3 * k * (3 * (c2[i] - cs) ** 2 * alpha[i]),  # c_j

                        -3 * k * (c2[i] - cs) ** 3,  # alpha_j
                        
                        0.5
                        * (D[i + 1] + D[i])
                        / (h2s * x2d_nodes[i] * x2dph),  # c_{j+1}
                        
                    ]
                    # residuals for Crank-Nicolson method
                    b += [
                        2 * (c2[i] - c2_o[i]) / dt
                        - 0.5
                        * (D[i + 1] + D[i])
                        * (c2[i + 1] - c2[i])
                        / (h2s * x2d_nodes[i] * x2dph)
                        + 0.5
                        * (D[i] + D[i - 1])
                        * (c2[i] - c2[i - 1])
                        / (h2s * x2d_nodes[i] * x2dmh)
                        + 3 * k * ((c2[i] - cs) ** 3) * alpha[i]
                        # old values for Crank-Nicolson
                        - 0.5
                        * (D_o[i + 1] + D_o[i])
                        * (c2_o[i + 1] - c2_o[i])
                        / (h2s * x2d_nodes[i] * x2dph)
                        + 0.5
                        * (D_o[i] + D_o[i - 1])
                        * (c2_o[i] - c2_o[i - 1])
                        / (h2s * x2d_nodes[i] * x2dmh)
                        + 3 * k * ((c2_o[i] - cs) ** 3) * alpha_o[i]
                    ]

    #
                ########################
                # CONSUMPTION OF U EQN #
                ########################
                #
            if c2[i] > cs:
                # [alpha] equation, reaction is active
                #
                #
                row += 2 * [j + 1]
                col += [j, j + 1]  # c_{j},alpha_j
                val += [
                    -eps * k * (3 * (c2[i] - cs) ** 2 * alpha[i]),
                    -eps * k * (c2[i] - cs) ** 3 - 2 / dt,
                ]
                # residuals for a Crank-Nicolson method
                b += [
                    2 * (alpha[i] - alpha_o[i]) / dt
                    + eps * k * ((c2[i] - cs) ** 3) * alpha[i]
                    + eps * k * ((c2_o[i] - cs) ** 3) * alpha_o[i]
                ]
            else:
                # [alpha] equation, no reaction terms
                #
                #
                row += [j + 1]
                col += [j + 1]  # alpha_j
                val += [-1 / dt]
                b += [(alpha[i] - alpha_o[i]) / dt]
        


        # at the last node we impose no H flux out of the domain
        i = n2 - 1
        j = n3 + ej * i
        row += 3 * [j]
        col += [j - 2 * ej, j - ej, j]
        val += [-1, 4, -3]
        # second order differencing
        b += [(3 * c2[i] - 4 * c2[i - 1] + c2[i - 2])]


        
        # solve the matrix problem for this iteration
        a = sp.sparse.coo_matrix((val, (row, col)), shape=(n3 + n2 * ej, n3 + n2 * ej))

        print(len(val),len(col),len(row),len(b))
        # system = petsc.PETScSparseLinearSystem(a, b)
        # x = system.linear_solve()

        # print('det',np.linalg.det(a.toarray()))
        
        x = sp.sparse.linalg.spsolve(a.tocsr(),b)
    
        # residual is the maximum correction value
        res = sp.linalg.norm(x, ord=np.inf)
        print(res)
        
        # print(
        #    f"iteration = {iteration} residual = {residual} ||b||_inf = {sp.linalg.norm(b, ord=np.inf)}"
        # )
        # add the corrections to the current guess
        c3[0:n3] += x[0:n3]
        c2[0:n2] += x[n3 : n3 + n2 * ej : ej]
        alpha[0:n2] += x[n3 + 1 : n3 + n2 * ej : ej]
        D = D1 * (1-alpha) + D2 * alpha

        if iterations > 10:
            print("Too many iterations")
            sys.exit()