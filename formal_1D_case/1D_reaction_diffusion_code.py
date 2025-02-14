import numpy as np
import scipy as sp

# Crank-Nicholson approach to solving reaction diffusion equation in 1 dimension iterating

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