import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

plt.style.use('formal_work/minorticks.mplstyle')



# some functions to remap coords
def x_from_X(X):
    """Defines the (non-uniformly distributed) physical coordinate 'x'
    in terms of the (uniformly distributed) computational coordinate
    'X'. Returns: not only x(X) but also x'(X) and x''(X).

    """
    BX = 0.001  # smaller BX leads to more NON-uniformity
    # +0*X below prevents numpy from returning a scalar when acting
    # on an iterable container for X
    return (X + BX) * (X + BX) - BX * BX, 2 * (X + BX)


# indexing is 1=Hydride (not present), 2=Metal, 3=Oxide
#
# discretisation

n2 = 16001  # metal
# number of variables in the bulk
nv = 3  # [H] and [U] in the bulk plus diffusivity of mixture

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 5e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
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

Drefstar = D2star  # using the U value

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

kstar = 2200
m = 1   # reaction order

reactK = (
    Lrefstar**2*N2star**m * eps**(m-1)/Drefstar * kstar  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
)


# new_scalings




# tscaling = (1/(eps*reactK))
# Lrefstar = (D2star/reactK)**0.5

tscaling = Lrefstar**2/D2star


print('x scaling : ', Lrefstar)
print('t scaling : ', tscaling)



D = D2 * np.ones(n2, dtype=np.float64)  # mixture diffusion in the bulk, initially =D2

# information output sanity check
print(f"Diffusivity ratio Dref*={D1star/Drefstar}")

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
dt = 0.000001
tol = 1.0e-8
threshold_alpha = 0.97


print(x2_nodes[200])

# time step the system
t = 0
step = 1  # periodic output counter

hydride_interface097 = []
hydride_interface098 = []
hydride_interface099 = []

V_list = []
t_list = []

numerics_lenth_hydride = []
asymptotics_length_hydride = []


dt_spacing = dt * 50

t_values = [dt_spacing*i+0.0000001 for i in range(1,int(150))]
t_values = np.array(t_values)


for t in t_values:
    alpha_tot = np.loadtxt(
            f'formal_work/data1D/c2_condon_097{t:.6f}.dat',
            )

    alpha = alpha_tot[:,1]
    x2_nodes = alpha_tot[:,0]

    interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - 0.97)
    roots = interpolated_alpha.roots(extrapolate=False)

    if len(roots) > 0:
        hydride_interface097.append(roots[0])
    else:
        hydride_interface097.append(0)

for t in t_values:
    alpha_tot = np.loadtxt(
            f'formal_work/data1D/c2_condon_098{t:.6f}.dat',
            )

    alpha = alpha_tot[:,1]
    x2_nodes = alpha_tot[:,0]

    interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - 0.98)
    roots = interpolated_alpha.roots(extrapolate=False)

    if len(roots) > 0:
        hydride_interface098.append(roots[0])
    else:
        hydride_interface098.append(0)

for t in t_values:
    alpha_tot = np.loadtxt(
            f'formal_work/data1D/c2_condon_099{t:.6f}.dat',
            )

    alpha = alpha_tot[:,1]
    x2_nodes = alpha_tot[:,0]

    interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - 0.99)
    roots = interpolated_alpha.roots(extrapolate=False)

    if len(roots) > 0:
        hydride_interface099.append(roots[0])
    else:
        hydride_interface099.append(0)
    

asymptotics_length_hydride = []
asymptotics_length_hydride.append(np.sqrt(2*D2star*kstar*(Castar-Csstar)**(m+1)/(3*N2star*(m+1)))/(1-0.97)*t_values * Lrefstar**2/Drefstar)
asymptotics_length_hydride.append(np.sqrt(2*D2star*kstar*(Castar-Csstar)**(m+1)/(3*N2star*(m+1)))/(1-0.98)*t_values * Lrefstar**2/Drefstar)
asymptotics_length_hydride.append(np.sqrt(2*D2star*kstar*(Castar-Csstar)**(m+1)/(3*N2star*(m+1)))/(1-0.99)*t_values * Lrefstar**2/Drefstar)




# print(asymptotics_length_hydride)

plt.plot(t_values*Lrefstar**2/Drefstar,np.array(hydride_interface097)*Lrefstar*1e7,color='red',label=r'$\alpha_c=0.97$',linestyle='dashed')
plt.plot(t_values*Lrefstar**2/Drefstar,asymptotics_length_hydride[0] * 1e7,color='red',)

plt.plot(t_values*Lrefstar**2/Drefstar,np.array(hydride_interface098)*Lrefstar*1e7,color='blue',label=r'$\alpha_c=0.98$',linestyle='dashed')
plt.plot(t_values*Lrefstar**2/Drefstar,asymptotics_length_hydride[1] * 1e7,color='blue')

plt.plot(t_values*Lrefstar**2/Drefstar,np.array(hydride_interface099)*Lrefstar*1e7,color='green',label=r'$\alpha_c=0.99$',linestyle='dashed')
plt.plot(t_values*Lrefstar**2/Drefstar,asymptotics_length_hydride[2] * 1e7,color='green')




# plt.xlabel(r'$\text{time / }\mathrm{s}$')
# plt.ylabel(r'$\Gamma^* \text{ / }\mathrm{nm}$')

plt.savefig('formal_work/1D_code/regime2a_plots/figures_latex/regime3-V-plot.eps',format='eps',bbox_inches='tight')

# plt.legend()
plt.show()









