import numpy as np
import scipy as sp

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
n2 = 101
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
alpha = np.zeros(n2)

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

k = 1e6

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
    c2_o = np.copy(c2)
    c3_o = np.copy(c3)
    alpha_o = np.copy(alpha)

    res = 1
    while res > tolerance:
        # Computes diffusivities
        D = D1 * (1-alpha) + D2 * alpha
        D_o = D1 * (1-alpha_o) + D2 * alpha_o
        row = []
        col = []
        val = []
        b = []

        # gas-oxide interface
        row += [0]
        col += [0]
        val += [1]
        b += 1 - c3[0] # c + c_tilde = 1 at surface

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
        
        for i in range(1,n2-1):
            j = n3 + ej * j # scales up the number of new matrix columns/rows noting that we need extra entries for alpha terms
            D = D1 * (1-alpha) + D2 * alpha

            if c2[i] < cs:  # no reaction
                x2dph = 0.5 * (x2d_nodes[j + 1] + x2d_nodes[j])
                x2dmh = 0.5 * (x2d_nodes[j] + x2d_nodes[j - 1])
                




        # at the last node we impose no H flux out of the domain
        j = n2 - 1
        k = n3 + ej * j
        row += 3 * [k]
        col += [k - 2 * ej, k - ej, k]
        val += [-1, 4, -3]
        # second order differencing
        b += [(3 * c2[j] - 4 * c2[j - 1] + c2[j - 2])]