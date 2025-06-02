import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
t_values = [1000*i+0.00001 for i in range(1,int(21))]
t_values = np.array(t_values)

t_05 = np.array([500*i+0.0001 for i in range(1,int(30))])

L2star = 500 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 20 * 1.0e-7  # oxide domain lengthscale 10nm
Lmstar = 500 * 1.0e-7 # 100um
D1star = 1.0e-15  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 5.0e-10  # cm^2/s, diffusion of H in U (room temp value)
D3star = 1.0e-13  # cm^2/s, diffusion of H in UO2 (room temp value)
N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 1.0e-5  # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castar = 1.0e-4  # mol/cm^3, surface value for [H] -- *ad-hoc* value

# fixed reference values to keep time/lengthscales consistent in comparisons
Lrefstar = 1.0e2 * 1.0e-7  # using 1um=1000nm here
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

r = 15 * 1e3 * 1e-7 / Lrefstar # radius

def ierfc(x):
    return 1/np.sqrt(np.pi) * np.exp(-x**2) - x * (1 - sp.special.erf(x))
    
t_star = 1000

non_dim_t = t_star * Drefstar / Lrefstar**2

t_dim = (t_values * Lrefstar**2 / Drefstar)

t_dim05 = (t_05 * Lrefstar**2 / Drefstar)


def c(t,z):
    return 2 * D3 * (t/D2)**0.5 / L3 * (ierfc(z/(4*D2*t)**0.5)-ierfc(((r**2 + z**2)/(4*D2*t))**0.5))

print(c(10000,0))

glascott0 = []
glascott5 = []
glascott10 = []
glascott15 = []

numerics0 = []
numerics5 = []
numerics10 = []
numerics15 = []

numerics01 = []
numerics51 = []
numerics101 = []
numerics151 = []

numerics0t05 = []
numerics5t05 = []
numerics10t05 = []
numerics15t05 = []

x2_nodes = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/x2_nodes_only_7.dat')

x2_nodes_l = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/x2_nodes_only_7.dat')
print(x2_nodes)
for t in t_values:
    c2 = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/c2_401_only_diff_7{t:.2f}.dat')
    c21 = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/c2_401_only_diff_8{t:.2f}.dat')
    c2_0 = c2[:,0]
    c2_01 = c21[:,0]
    

    spline = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,c2_0,k=4)

    spline1 = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes_l,c2_01,k=4)


    glascott0.append(c(t,0))
    glascott5.append(c(t,5*1e3*1e-7/Lrefstar))
    glascott10.append(c(t,10*1e3*1e-7/Lrefstar))
    glascott15.append(c(t,15*1e3*1e-7/Lrefstar))
    
    numerics0.append(spline(0))
    numerics5.append(spline(5*1e3*1e-7/Lrefstar))
    numerics10.append(spline(10*1e3*1e-7/Lrefstar))
    numerics15.append(spline(15*1e3*1e-7/Lrefstar))

    numerics01.append(spline1(0))
    numerics51.append(spline1(5*1e3*1e-7/Lrefstar))
    numerics101.append(spline1(10*1e3*1e-7/Lrefstar))
    numerics151.append(spline1(15*1e3*1e-7/Lrefstar))


# for t in t_05:
#     c2_t05 = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/c2_201_t05{t:.2f}.dat')
    
#     c2_t05 = c2_t05[:,0]
    
#     splinet05 = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,c2_t05,k=4)

#     numerics0t05.append(splinet05(0))
#     numerics5t05.append(splinet05(5*1e3*1e-7/Lrefstar))
#     numerics10t05.append(splinet05(10*1e3*1e-7/Lrefstar))
#     numerics15t05.append(splinet05(15*1e3*1e-7/Lrefstar))



plt.plot(t_dim,glascott0,color='red',label='z=0')
plt.plot(t_dim,glascott5,color='blue',label='z=5um')
plt.plot(t_dim,glascott10,color='orange',label='z=10um')
plt.plot(t_dim,glascott15,color='green',label='z=15um')

plt.plot(t_dim,numerics0,color='red',linestyle='dashed')
plt.plot(t_dim,numerics5,color='blue',linestyle='dashed')
plt.plot(t_dim,numerics10,color='orange',linestyle='dashed')
plt.plot(t_dim,numerics15,color='green',linestyle='dashed')

plt.plot(t_dim,numerics01,color='red',linestyle='dashdot')
plt.plot(t_dim,numerics51,color='blue',linestyle='dashdot')
plt.plot(t_dim,numerics101,color='orange',linestyle='dashdot')
plt.plot(t_dim,numerics151,color='green',linestyle='dashdot')


# plt.plot(t_dim05,numerics0t05,color='red',linestyle='dotted')
# plt.plot(t_dim05,numerics5t05,color='blue',linestyle='dotted')
# plt.plot(t_dim05,numerics10t05,color='orange',linestyle='dotted')
# plt.plot(t_dim05,numerics15t05,color='green',linestyle='dotted')




plt.legend()
plt.show()