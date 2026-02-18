
"""e.g. python diffusionBC_reaction.py -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist"""

import numpy as np
import scipy as sp
# import deps.sp2petsc as petsc
import sys
from matplotlib import pyplot as plt
import time 

plt.style.use('formal_work/minorticks.mplstyle')

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

n2 = 1001  # metal
# number of variables in the bulk
nv = 3  # [H] and [U] in the bulk plus diffusivity of mixture


tempC = 75
pressure = 1


tempK = 273.15 + tempC

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 1e2 * 1e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 7.5 * 1.0e-7  # oxide domain lengthscale 10nm
D1star = 1.0e-13  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 1.94e-2 * np.exp(-5570/tempK)

D3star = 0.037 * np.exp(-7200/tempK) # cm^2/s, diffusion of H in UO2

N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 0.2446 * np.exp(-52000 / (8.31 * tempK)) # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castar = N2star * 4.13e-6 * np.exp(-894/tempK) * np.sqrt(pressure)  # mol/cm^3, surface value for [H] -- *ad-hoc* value

kstar = 10.4/N2star * np.exp(1542/tempK) # reaction constant in cm^3/mol


# fixed reference values to keep time/lengthscales consistent in comparisons
# Lrefstar = 100 * 1e-7  # using 1um=1000nm here
Lrefstar = L2star



D_factor = 10000  # Default value is 1e4

D_largestar = D_factor * D2star
Drefstar = D_largestar # using the U value

D_large = D_largestar / Drefstar

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



m = 1   # reaction order

reactK = (
    1  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
)

Lrefstar = np.sqrt( Drefstar / (N2star * kstar))


# tscaling = (1/(eps*reactK))
# Lrefstar = (D2star/reactK)**0.5

# lambda = 0.1 * pressure

tscaling = Lrefstar**2/Drefstar / eps


print('x scaling : ', Lrefstar)
print('t scaling : ', tscaling)


values = ['1','05','025']
legend_values = ['press=1','press=0.5','press=0.25']
multiplier = [1,0.5,0.25]

legend = dict(zip(values, legend_values))
multiplier_dic = dict(zip(values, multiplier))
dic_alpha_int = {el:[] for el in values}

print(dic_alpha_int)


for i in values:
    t_values = np.loadtxt(
        f'formal_work/data1D/lambdaP/tvalues_D1e4_pressure{i}_75C.dat',
        )
    for t in t_values:
        c2_tot = np.loadtxt(
        f'formal_work/data1D/lambdaP/c2_condon_oxide_D1e4_pressure{i}_75C_{t:.7f}.dat',
        )
        alpha_tot = np.loadtxt(
            f'formal_work/data1D/lambdaP/alpha_condon_oxide_D1e4_pressure{i}_75C_{t:.7f}.dat',
            )
        
        x2_nodes = c2_tot[:,0]
        c2 = c2_tot[:,1]
        alpha = alpha_tot[:,1]
    
        dic_alpha_int[i].append(alpha[0])

    plt.plot(t_values*tscaling,dic_alpha_int[i],label=legend[i])
    plt.ylabel(r'$\alpha$ at interface')
    plt.xlabel('time (seconds)')
plt.legend()
plt.show()