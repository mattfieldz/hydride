import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
t_values = [0.01*i+0.00001 for i in range(1,int(20))]
t_values = np.array(t_values)

t_05 = np.array([500*i+0.0001 for i in range(1,int(30))])

L2star = 1 * 1e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 20 * 1.0e-7  # oxide domain lengthscale 10nm
Lmstar = 150 * 1e3 * 1.0e-7 # 100um
D1star = 1.0e-13  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 1.49e-10  # cm^2/s, diffusion of H in U (room temp value)
D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 (room temp value)
N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 1.0e-5  # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castar = 1.0e-4  # mol/cm^3, surface value for [H] -- *ad-hoc* value

# fixed reference values to keep time/lengthscales consistent in comparisons
Lrefstar = 100 * 1e-7  # using 1um=1000nm here
Lrefstar = L2star
Drefstar = D2star  # using the U value

# non-dimensional domains
L3 = L3star / Lrefstar  # oxide
L2 = L2star / Lrefstar  # bulk
Lm = Lmstar / Lrefstar # width

D3 = D3star / Drefstar
D2 = D2star / Drefstar
D1 = D1star / Drefstar
cs = Csstar / Castar  # saturation value \in [0,1)
eps = Castar / N2star  # relative forcing measure
reactK = (
    0  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
)

r = 1e3 * 1e-7 / Lrefstar # radius

def ierfc(x):
    return 1/np.sqrt(np.pi) * np.exp(-x**2) - x * (1 - sp.special.erf(x))
    
t_star = 1000

non_dim_t = t_star * Drefstar / Lrefstar**2

t_dim = (t_values * Lrefstar**2 / Drefstar)

t_dim05 = (t_05 * Lrefstar**2 / Drefstar)


def c(t,z):
    return 2 * D3 * (t/D2)**0.5 / L3 * (ierfc(z/(4*D2*t)**0.5)-ierfc(((r**2 + z**2)/(4*D2*t))**0.5))

t_values_dim = np.linspace(0,10,1000)

t_values = 1/(Lrefstar ** 2 / Drefstar ) * t_values_dim

c2 = c(t_values,0)

print(t_dim,c2)