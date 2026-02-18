import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
t_values = [1000*i+0.00001 for i in range(1,int(50))]
t_values = np.array(t_values)

t_05 = np.array([30*i+0.000001 for i in range(1,int(50))])

L2star = 150 * 1e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 20 * 1.0e-7  # oxide domain lengthscale 10nm
Lmstar = 150 * 1e3 * 1.0e-7 # 100um
D1star = 1.0e-13  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 5.0e-10  # cm^2/s, diffusion of H in U (room temp value)
D3star = 1.0e-13  # cm^2/s, diffusion of H in UO2 (room temp value)
N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 1.0e-6  # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castar = 1.0e-4  # mol/cm^3, surface value for [H] -- *ad-hoc* value

# fixed reference values to keep time/lengthscales consistent in comparisons
Lrefstar = 100 * 1e-7
  # using 1um=1000nm here
# Lrefstar = L2star
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

numerics0_16 = []
numerics5_16 = []
numerics10_16 = []
numerics15_16 = []

numerics0_asy = []
numerics5_asy = []
numerics10_asy = []
numerics15_asy = []

numerics0t05 = []
numerics5t05 = []
numerics10t05 = []
numerics15t05 = []

k14 = Lrefstar**2*eps**2*N2star**3/Drefstar * 1e14
k16 = Lrefstar**2*eps**2*N2star**3/Drefstar * 1e16
k12 = Lrefstar**2*eps**2*N2star**3/Drefstar * 1e12
k14_05 = L2star**2*eps**2*N2star**3/Drefstar * 1e14
k16_05 = L2star**2*eps**2*N2star**3/Drefstar * 1e16
k17 = Lrefstar**2*eps**2*N2star**3/Drefstar * 1e17
print(k17)
print(k16)


t_dim05 = (t_05 * (Lrefstar)**2 / (Drefstar))
# print(t_dim05)
# print(t_dim)

x2_nodes = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/x2_nodes.dat')

x2_nodes_l = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/x2_nodes.dat')
# print(x2_nodes)
t_c = 0
for t in t_values:
    
    c2_14_tot = np.loadtxt(f'formal_work/data1D/k1e17/c2{t:.2f}.dat')
    c2_16_tot = np.loadtxt(f'formal_work/data1D/k1e16/c2{t:.2f}.dat')
    
    a_14_tot = np.loadtxt(f'formal_work/data1D/k1e17/alpha{t:.2f}.dat')
    a_16_tot = np.loadtxt(f'formal_work/data1D/k1e16/alpha{t:.2f}.dat')
    
    alpha_asy_tot = np.loadtxt(f'formal_work/data1D/asymptotic/alpha_1001_{t_05[t_c]:.5f}.dat')
    c2_asy_tot = np.loadtxt(f'formal_work/data1D/asymptotic/c2_1001_{t_05[t_c]:.5f}.dat')

    x2_nodes_asy = c2_asy_tot[:,0]
    # c2_asy = cs + c2_asy_tot[:,1] * k12**(-0.25)
    c2_asy = c2_asy_tot[:,1] 
    a_asy = alpha_asy_tot[:,1]
    x2_nodes_14 = c2_14_tot[:,0]

    c2_14 = (c2_14_tot[:,1] - cs) * k17**0.25
    c2_16 = (c2_16_tot[:,1] - cs) * k16**0.25

    a_14 = (a_14_tot[:,1])
    a_16 = (a_16_tot[:,1]) 
    # print(a_14)
    # c2_14 = c2_14_tot[:,1]
    # c2_16 = c2_16_tot[:,1]

    x2_nodes_16 = c2_16_tot[:,0]
    # c2_01 = c21[:,0]
    # print(c2_0)
    
    spline_14 = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes_14,a_14,k=4)

    spline_asy = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes_asy,a_asy,k=4)

    spline_16 = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes_16,a_16,k=4)

    glascott0.append(c(t,0))
    glascott5.append(c(t,5*1e3*1e-7/Lrefstar))
    glascott10.append(c(t,10*1e3*1e-7/Lrefstar))
    glascott15.append(c(t,15*1e3*1e-7/Lrefstar))
    
    numerics0.append(spline_14(0))
    numerics5.append(spline_14(5*1e3*1e-7/Lrefstar))
    numerics10.append(spline_14(10*1e3*1e-7/Lrefstar))
    numerics15.append(spline_14(15*1e3*1e-7/Lrefstar))

    numerics0_asy.append(spline_asy(0))
    numerics5_asy.append(spline_asy(5*1e3*1e-7/L2star))
    numerics10_asy.append(spline_asy(10*1e3*1e-7/L2star))
    numerics15_asy.append(spline_asy(15*1e3*1e-7/L2star))

    numerics0_16.append(spline_16(0))
    numerics5_16.append(spline_16(5*1e3*1e-7/Lrefstar))
    numerics10_16.append(spline_16(10*1e3*1e-7/Lrefstar))
    numerics15_16.append(spline_16(15*1e3*1e-7/Lrefstar))

    t_c += 1
# for t in t_05:
#     c2_t05 = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/c2_201_t05{t:.2f}.dat')
    
#     c2_t05 = c2_t05[:,0]
    
#     splinet05 = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,c2_t05,k=4)

#     numerics0t05.append(splinet05(0))
#     numerics5t05.append(splinet05(5*1e3*1e-7/Lrefstar))
#     numerics10t05.append(splinet05(10*1e3*1e-7/Lrefstar))
#     numerics15t05.append(splinet05(15*1e3*1e-7/Lrefstar))



# plt.plot(t_dim,glascott0,color='red',label='z=0')
# plt.plot(t_dim,glascott5,color='blue',label='z=5um')
# plt.plot(t_dim,glascott10,color='orange',label='z=10um')
# plt.plot(t_dim,glascott15,color='green',label='z=15um')

plt.plot(t_dim*k17**0.25*eps,numerics0,color='red',linestyle='dashed')
# plt.plot(t_dim,numerics5,color='blue',linestyle='dashed')
# plt.plot(t_dim,numerics10,color='orange',linestyle='dashed')
# plt.plot(t_dim,numerics15,color='green',linestyle='dashed')

plt.plot(t_dim05,numerics0_asy,color='red',linestyle='solid')

# plt.plot(t_dim05,numerics5_asy,color='blue',linestyle='dashdot')
# plt.plot(t_dim05,numerics10_asy,color='orange',linestyle='dashdot')
# plt.plot(t_dim05,numerics15_asy,color='green',linestyle='dashdot')

# plt.plot(t_dim*k16**0.25*eps,numerics0_16,color='red',linestyle='dotted')
# plt.plot(t_dim,numerics5_16,color='blue',linestyle='dotted')
# plt.plot(t_dim,numerics10_16,color='orange',linestyle='dotted')
# plt.plot(t_dim,numerics15_16,color='green',linestyle='dotted')

# plt.plot(t_dim05,numerics0t05,color='red',linestyle='dotted')
# plt.plot(t_dim05,numerics5t05,color='blue',linestyle='dotted')
# plt.plot(t_dim05,numerics10t05,color='orange',linestyle='dotted')
# plt.plot(t_dim05,numerics15t05,color='green',linestyle='dotted')

plt.ylabel('$c$ / $C^*$')
plt.xlabel('t / s')

plt.legend()
plt.show()