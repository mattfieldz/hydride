from matplotlib import pyplot as plt
import numpy as np
import scipy as sp


plt.style.use('formal_work/minorticks.mplstyle')

dt = 0.03

# t_values = [dt*i+0.00001 for i in range(1,int(29.40/dt))]
t_values = [dt*i+0.00001 for i in range(1,int(13.86/dt))]

t_values = np.array(t_values)

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 1e2 * 1e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
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

D_factor = 1e4

Drefstar = D_factor * D2star  # using the U value

# Drefstar = D2star

D_large = 1 

D1 = D1star/Drefstar
D2 = D2star/Drefstar
D3 = D3star/Drefstar

# non-dimensional domains
L3 = L3star / Lrefstar  # oxide
L2 = L2star / Lrefstar  # bulk

cs = Csstar / Castar
eps = Castar / N2star

print(eps,'eps')

kstar = 1e4
m = 1   # reaction order

reactK = (
    Lrefstar**2*N2star/Drefstar * kstar  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
)

print('react_k',reactK)



threshold_alpha = 0.98

t_scale_for_alpha_c = eps * (reactK / 2 * (1-threshold_alpha**2) + reactK**1.5 * (2*L3 / (5*D3)) * (3*D2)**0.5 * (1-threshold_alpha**2.5))

print('lol',eps*reactK / 2 * (1-threshold_alpha**2))


n2 = 4001
x1_asy = []

hydride_interface = [0]
hydride_interface_D = [0]
step = 0

tscaling = Lrefstar**2 / Drefstar

print(tscaling,'tscaling')

t_values_dim = t_values * tscaling

print(t_scale_for_alpha_c,'non-dim')
print(t_scale_for_alpha_c * tscaling,'dim')


V = [0]
V_deriv = [0]
c2_int = []

V_D = [0]
V_deriv_D = [0]
c2_int_D = []

alpha_int = []
alpha_int_D = []

for t in t_values:
    
    c2_tot = np.loadtxt(
        f'formal_work/data1D/c2_condon_oxide_larger_D{t:.2f}.dat',
        )
    alpha_tot = np.loadtxt(
        f'formal_work/data1D/alpha_condon_oxide_larger_D{t:.2f}.dat',
        )
    
    x2_nodes = c2_tot[:,0]
    c2 = c2_tot[:,1]
    alpha = alpha_tot[:,1]
    
    interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - threshold_alpha)
    interpolated_c2 = sp.interpolate.CubicSpline(x2_nodes, c2)
    roots = interpolated_alpha.roots(extrapolate=False)
    if len(roots) > 0:
        hydride_interface.append(roots[0])
        
        
    else:
        hydride_interface.append(0)
    if step > 0:
        V.append(Lrefstar*((hydride_interface[step])-hydride_interface[step-1])/(dt*tscaling))

        V_deriv.append((V[step]-V[step-1])/(dt*tscaling))
     

    c2_int.append(c2[0])

    c2_tot_D = np.loadtxt(f"formal_work/data1D/c2_condon_oxide_Dfac1e3{0.1*t:.3f}.dat")  
    alpha_tot_D = np.loadtxt(f"formal_work/data1D/alpha_condon_oxide_Dfac1e3{0.1*t:.3f}.dat")


    x2_nodes = c2_tot_D[:,0]
    c2_D = c2_tot_D[:,1]
    alpha_D = alpha_tot_D[:,1]
    
    interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha_D - threshold_alpha)
    interpolated_c2 = sp.interpolate.CubicSpline(x2_nodes, c2_D)
    roots = interpolated_alpha.roots(extrapolate=False)
    if len(roots) > 0:
        hydride_interface_D.append(roots[0])
        
        
    else:
        hydride_interface_D.append(0)
    if step > 0:
        V_D.append(Lrefstar*((hydride_interface_D[step])-hydride_interface_D[step-1])/(dt*tscaling))

        V_deriv_D.append((V_D[step]-V_D[step-1])/(dt*tscaling))
    step += 1    

    c2_int_D.append(c2_D[0])

    alpha_int.append(alpha[0])
    alpha_int_D.append(alpha_D[0])
V_arr = np.array(V)
V_deriv_arr = np.array(V_deriv)
V_arr_D = np.array(V_D)
V_deriv_arr_D = np.array(V_deriv_D)
hydride_interface_arr = np.array(hydride_interface)
hydride_interface_arr_D = np.array(hydride_interface_D)

# plt.plot(t_values_dim[100:-1],1e7 * V_arr[100:-1],label='D=1e4')
# plt.plot(t_values_dim[100:-1],1e7 * V_arr_D[100:-1],label='D=1e3')

# plt.plot(t_values_dim[0:-1],1e7 * Lrefstar * hydride_interface_arr[1:-1],label='D=1e4')
# plt.plot(t_values_dim[0:-1],1e7 * Lrefstar * hydride_interface_arr_D[1:-1],label='D=1e3')
# plt.plot(t_values_dim[0:-1],1e7 * Lrefstar * np.sqrt(2 * 1e4 * eps * (1-cs) / (3) * (D2star / Lrefstar**2) * t_values_dim[0:-1]),)

# plt.plot(t_values_dim[0:100], np.array(c2_int)[0:100],label='D=1e4')
# plt.plot(t_values_dim[0:100], np.array(c2_int_D)[0:100],label='D=1e3')


# plt.plot(1e7 * Lrefstar * x2_nodes[0:100],c2[0:100],label='D=1e4')
# plt.plot(1e7 * Lrefstar * x2_nodes[0:100],c2_D[0:100],label='D=1e3')
# plt.plot(x2_nodes,c2_D)

D_arr = np.zeros(100)
D_arr_D = np.zeros(100)

for j in range(100): 
    if alpha[j] > threshold_alpha:
        D_arr[j] = 1
    elif alpha[j] > 0.8:
        D_arr[j] = 1 + 1e4 / (0.8-threshold_alpha)**3 * (alpha[j]-threshold_alpha)**3
    else:
        D_arr[j] = 1e4
        
    if alpha_D[j] > threshold_alpha:
        D_arr_D[j] = 1
    elif alpha_D[j] > 0.8:
        D_arr_D[j] = 1 + 1e3 / (0.8-threshold_alpha)**3 * (alpha_D[j]-threshold_alpha)**3
    else:
        D_arr_D[j] = 1e3

# plt.plot(1e7 * Lrefstar * x2_nodes[0:100],D_arr,label='D=1e4')
# plt.plot(1e7 * Lrefstar * x2_nodes[0:100],D_arr_D,label='D=1e3')

# plt.plot(t_values[250:-1] * tscaling,1e7 * V_deriv_arr[250:-1])
# plt.plot(t_values_dim[10:-1],c2_int[10:-1])
# plt.plot(t_values_dim[10:-1],c2_int_D[10:-1])
# plt.plot(x2_nodes,c2)
# plt.plot(Lrefstar * 1e7 * x2_nodes[0:100],alpha[0:100])

# plt.xlabel('z / nm')
# plt.ylabel(r'$D / D_2$')


# plt.legend()
# plt.show()

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.7)

twin1 = ax.twinx()



p1, = ax.plot(1e7 * Lrefstar * x2_nodes[0:100], alpha[0:100], "C0", label=r"$\alpha$")
p2, = twin1.plot(1e7 * Lrefstar * x2_nodes[0:100], c2[0:100], "C1", label=r"$c$")

ax.set(xlim=(0, 100),ylim=(0.5, 1),xlabel=r"$x^* / \mathrm{nm}$", ylabel=r"$\alpha$")
twin1.set(ylabel=r"$c$")

twin1.grid(False)

ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())


ax.tick_params(axis='y', colors=p1.get_color())
twin1.tick_params(axis='y', colors=p2.get_color())

print(t_values_dim[-1])

ax.legend(handles=[p1, p2])
plt.show()


   