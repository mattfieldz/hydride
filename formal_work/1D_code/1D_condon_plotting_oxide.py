from matplotlib import pyplot as plt
import numpy as np
import scipy as sp

def x_from_X(X):
    """Defines the (non-uniformly distributed) physical coordinate 'x'
    in terms of the (uniformly distributed) computational coordinate
    'X'. Returns: not only x(X) but also x'(X) and x''(X).

    """
    BX = 0.01  # smaller BX leads to more NON-uniformity
    # +0*X below prevents numpy from returning a scalar when acting
    # on an iterable container for X
    return (X + BX) * (X + BX) - BX * BX, 2 * (X + BX)




dt = 0.003

# t_values = [dt*i+0.00001 for i in range(1,int(29.40/dt))]
t_values = [dt*i+0.00001 for i in range(1,int(0.2880/dt))]

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

D_factor = 1e2

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

integral_quantity = [0,0]

c2_end = []
c2_end_D = []
c = 0

Xmax = sp.optimize.root(lambda X: x_from_X(X)[0] - L2, 5).x[0]
# exit if sanity check fails
assert abs(x_from_X(Xmax)[0] - L2) < 1.0e-8
# nodes in the oxide and bulk

X2_nodes = np.linspace(0, Xmax, n2)  # bulk nodes in comp. coordinate
# physical node locations
x2_nodes, x2d_nodes = x_from_X(X2_nodes)





for t in t_values[0:-1]:
    
    if c != 0:
        c2_prev = np.copy(c2)
        alpha_prev = np.copy(c2)
    c += 1
    c2_tot = np.loadtxt(
        f'formal_work/data1D/c2_condon_oxide_no_threshold{10*t:.3f}.dat',
        )
    alpha_tot = np.loadtxt(
        f'formal_work/data1D/alpha_condon_oxide_no_threshold{10*t:.3f}.dat',
        )
    
    x2_nodes = c2_tot[:,0]
    c2 = c2_tot[:,1]
    alpha = alpha_tot[:,1]
    
    interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - threshold_alpha)
    interpolated_c2 = sp.interpolate.CubicSpline(x2_nodes, c2-cs)

    just_alpha = sp.interpolate.CubicSpline(x2_nodes,alpha)


    roots = interpolated_c2.roots(extrapolate=False)
    if len(roots) > 0:
        hydride_interface.append(roots[0])
        
        integral = just_alpha.integrate(0,roots[0])
        integral_quantity.append(integral)
    else:
        hydride_interface.append(0)
    if step > 0:
        V.append(Lrefstar*((hydride_interface[step])-hydride_interface[step-1])/(dt*tscaling))

        V_deriv.append((V[step]-V[step-1])/(dt*tscaling))
     

    c2_int.append(c2[0])

    c2_tot_D = np.loadtxt(f"formal_work/data1D/c2_condon_oxide_no_threshold_D1e2{t:.4f}.dat")  
    alpha_tot_D = np.loadtxt(f"formal_work/data1D/alpha_condon_oxide_no_threshold_D1e2{t:.4f}.dat")


    x2_nodes = c2_tot_D[:,0]
    c2_D = c2_tot_D[:,1]
    alpha_D = alpha_tot_D[:,1]
    
    interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha_D - threshold_alpha)

    
    interpolated_c2_D = sp.interpolate.CubicSpline(x2_nodes, c2_D-cs)
    roots = interpolated_c2_D.roots(extrapolate=False)
    if len(roots) > 0:
        hydride_interface_D.append(roots[0])
        
        
    else:
        hydride_interface_D.append(0)

    
    interpolated_c2_end_D = sp.interpolate.CubicSpline(x2_nodes, c2_D-0.0001)
    roots = interpolated_c2_end_D.roots(extrapolate=False)
    if len(roots) > 0:
        c2_end_D.append(roots[0])
        
        
    else:
        c2_end_D.append(0)

    interpolated_c2_end = sp.interpolate.CubicSpline(x2_nodes, c2-0.0001)
    roots = interpolated_c2_end.roots(extrapolate=False)
    if len(roots) > 0:
        c2_end.append(roots[0])
        
        
    else:
        c2_end.append(0)

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



# B = L3 * np.sqrt(3*D2)/ D3 * reactK**(-0.5)
# glascott_time_dim = np.pi / 4 * D2star * (cs * L3star)**2 / (D3star)**2

# alpha_predicted =  -( eps * (t_values[0:300]-glascott_time_dim / tscaling) - B ) / (np.exp(1)*reactK**(-1) + B)
# alpha_predicted =  1-eps/(np.exp(1)*reactK**(-1) + B) * (t_values[0:300] - glascott_time_dim / (tscaling))

# print(alpha_predicted[0:30],t_values[0:30])

# alpha_decay = threshold_alpha * np.exp(eps *np.sqrt(3) * reactK * (np.array(c2_int_D[0:500])-cs) * (0.35-t_values[0:500]))

# plt.plot(t_values_dim[0:500], np.array(alpha_int_D)[0:500],label='Numerics solution')
# plt.plot(t_values_dim[0:50],alpha_predicted[0:50],label='Linear asymptotic prediction')
# plt.plot(t_values_dim[5:500],alpha_decay[5:500],label='Exponential decaying prediction')
# plt.plot(t_values_dim[0:500],np.exp(-reactK * 3 * t_values[0:500]))
# plt.xlabel(r'$t^*/s$')
# plt.ylabel(r'$\alpha$')

# plt.plot(t_values_dim[5:500], np.array(c2_int_D)[5:500]-cs,label='D=1e3')


# plt.plot(1e7 * Lrefstar * x2_nodes[0:100],alpha[0:100],label='D=1e4')
# plt.plot(1e7 * Lrefstar * x2_nodes[0:100],alpha_D[0:100],label='D=1e3')
# plt.plot(x2_nodes,c2_D)

D_arr = np.zeros(302)
D_arr_D = np.zeros(302)

for j in range(302): 
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
# plt.plot(t_values_dim[50:-1],np.array(c2_int)[50:-1]-cs)
r_rate = (3*reactK*(c2-cs)*alpha)
for i in range(len(r_rate)):
    if r_rate[i] < 0:
        r_rate[i] = 0


diff_term = np.zeros(300)


h2 = L2 / (Xmax * (n2 - 1))
h2s = h2 * h2
print(reactK)
D_arr = D_arr / 1e4
D_arr_D = D_arr_D / 1e3
# D_arr = 1*np.ones(200)

c2_deriv = np.zeros(300)
c2_second_deriv = np.zeros(300)

# for i in range(300):
#     j = i + 1
#     x2dph = 0.5 * (x2d_nodes[j + 1] + x2d_nodes[j])
#     x2dmh = 0.5 * (x2d_nodes[j] + x2d_nodes[j - 1])

#     diff_term[i] = -0.5 * (D_arr_D[j + 1] + D_arr_D[j]) * (c2_D[j + 1] - c2_D[j]) / (h2s * x2d_nodes[j] * x2dph) 
#     + 0.5 * (D_arr[j] + D_arr[j - 1]) * (c2_D[j] - c2_D[j - 1]) / (h2s * x2d_nodes[j] * x2dmh)
#     # print(diff_term[i])

#     c2_deriv[i] = (c2_D[j+1] - c2_D[j-1]) / (2*h2 * x2d_nodes[j])
#     c2_second_deriv[i] = (c2_D[j+1] - c2_D[j]) / (h2s * x2dph * x2d_nodes[j]) + (c2_D[j-1] - c2_D[j]) / (h2s * x2dmh * x2d_nodes[j]) 



# plt.plot(x2_nodes[0:300] * Lrefstar * 1e7, (c2[0:300] - c2_prev[0:300])/(0.1*dt))
# plt.plot(x2_nodes[0:300] * Lrefstar * 1e7, r_rate[0:300])
# plt.plot(x2_nodes[0:300] * Lrefstar * 1e7, diff_term)
# plt.plot(x2_nodes[0:300] * Lrefstar * 1e7 ,diff_term/r_rate[0:300])


# plt.plot(x2_nodes[0:300],D_arr_D[0:300])
# plt.plot(x2_nodes[0:300], c2_second_deriv[0:300])
# plt.plot(x2_nodes[0:500], c2[0:500])

# plt.plot(t_values[250:-1] * tscaling,1e7 * V_deriv_arr[250:-1])
# plt.plot(t_values_dim[30:-1],np.array(c2_int[30:-1]))
# plt.plot(t_values_dim[30:-1],np.array(c2_int_D[30:-1]))

# plt.plot(t_values_dim[100:-1],(np.array(c2_int[100:-1])-np.array(c2_int[99:-2]))/(t_values_dim[1]-t_values_dim[0]))
tau_values = t_values * eps*reactK**0.5

print('log',D3 * (1-cs) * np.sqrt(np.pi/(D2*eps)) * reactK**(-0.25))

Gamma_approx = np.sqrt(4* D_large / eps) * np.sqrt(np.log( D3 * (1-cs) / (L3) * np.sqrt(np.pi/(D2*eps)) * reactK**(-0.25) )
                                                  *tau_values
                                                     + 0.5 * tau_values * (np.log(tau_values)-1) ) 

plt.plot(t_values_dim,1e7 * Lrefstar * hydride_interface_arr,label='Numerics,D=1e3')
plt.plot(t_values_dim,1e7 * Lrefstar * hydride_interface_arr_D,label='Numerics,D=1e2')
# plt.plot(t_values_dim[0:-1],alpha_int,label='numerics alpha D1e3')
# plt.plot(t_values_dim[0:-1],alpha_int_D,label='numerics alpha D1e2')

# plt.plot(x2_nodes[0:100],alpha[0:100])
# plt.plot(x2_nodes[0:100],alpha_D[0:100])


# plt.plot((t_values_dim*(np.log(t_values_dim)-1))**0.5,1e7 * Lrefstar * hydride_interface_arr,label='log behaviour?')

# plt.plot(t_values_dim,1e4*Lrefstar * Gamma_approx * reactK**(-0.5))


# plt.plot(t_values_dim,1e7 * Lrefstar * hydride_interface_arr_D)
# plt.plot(t_values_dim, 1e7 * Lrefstar * (np.array(integral_quantity) + eps / 3 * (D3/L3 * (1-cs) * t_values-2*cs/(np.sqrt(np.pi))*np.sqrt(D2*t_values))),label='Asymptotic approximation')
# plt.plot(t_values_dim, 1e7 * Lrefstar * (np.array(integral_quantity)))
# plt.plot(t_values[100:-1],np.array(c2_int[99:-1])-cs)
# plt.plot(t_values[100:-1],(np.array(c2_int[99:-1])-np.array(c2_int[98:-2]))/(t_values[4]-t_values[3]))

# plt.plot(t_values_dim,1e7 * Lrefstar * np.sqrt(2*eps*(1e-4)/(3*(1-threshold_alpha))*t_values))
# plt.plot(Lrefstar * 1e7 * x2_nodes[0:100],alpha[0:100])

plt.xlabel(r'$time^*/s$')
plt.ylabel(r'Hydride thickness / nm$')


plt.legend()
plt.show()

# fig, ax = plt.subplots()
# fig.subplots_adjust(right=0.7)

# twin1 = ax.twinx()


# p1, = ax.plot(1e7 * Lrefstar * x2_nodes[0:100], alpha_D[0:100], "C0", label=r"$\alpha$")
# p2, = twin1.plot(1e7 * Lrefstar * x2_nodes[0:100], c2_D[0:100], "C1", label=r"$c_2$")


# ax.set(xlim=(0, 100),ylim=(0.5, 1),xlabel=r"$z / \mathrm{nm}$", ylabel=r"$\alpha$")
# twin1.set(ylabel=r"$c_2$")


# ax.yaxis.label.set_color(p1.get_color())
# twin1.yaxis.label.set_color(p2.get_color())


# ax.tick_params(axis='y', colors=p1.get_color())
# twin1.tick_params(axis='y', colors=p2.get_color())

# print(t_values_dim[-1])

# ax.legend(handles=[p1, p2])
# plt.show()


   