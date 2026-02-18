from matplotlib import pyplot as plt
import numpy as np
import scipy as sp

plt.style.use('formal_work/minorticks.mplstyle')

dt = 20

t_values = [dt*i+0.00001 for i in range(1,int(50))]
t_values = np.array(t_values)

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 1e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
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

D1 = D1star/Drefstar
D2 = D2star/Drefstar
D3 = D3star/Drefstar

# non-dimensional domains
L3 = L3star / Lrefstar  # oxide
L2 = L2star / Lrefstar  # bulk

print(L3,'L3')

cs = Csstar / Castar
eps = Castar / N2star
back_stepping = np.loadtxt(f'1Dexamples/backtimestepsk7.dat')

Trefstar = Lrefstar ** 2 / Drefstar

k_values = [1e13,1e14,1e15,1e16]
k_s = ['k1e13','k1e14','k1e15','k1e16']

# k_values = [1e15,1e16]
# k_s = ['k1e15','k1e16']

x1 = [[] for i in range(len(k_values))]

t_induct = [0 for i in range(len(k_values))]
alpha_int = [[] for i in range(len(k_values))]

n2 = 1001
x1_asy = []

for i in range(len(k_values)):
    for t in t_values:
        
        c2_tot = np.loadtxt(f'formal_work/data1D/{k_s[i]}/c2{t:.2f}.dat')
        alpha_tot = np.loadtxt(f'formal_work/data1D/{k_s[i]}/alpha{t:.2f}.dat')
        x2_nodes = c2_tot[:,0]
        c2 = c2_tot[:,1]
        alpha = alpha_tot[:,1]

        # maximise reaction term
        spline = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,(c2-cs)**3*alpha,k=4)
        cr_points = spline.derivative().roots()

        # c=cs
        # spline = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,c2-cs,k=3)
        # cr_points = spline.roots()
        
        cr_values = spline(np.append(cr_points,(0,L2)))
        cr_maximal_index = np.argmax(cr_values)
        
        
        if len(cr_points) > 0:
            x1[i].append(np.append(cr_points,(0,L2))[cr_maximal_index])
            # x1[i].append(cr_points[0])
            if x1[i][-1] == 0:
                t_induct[i] += dt
        else:
            x1[i].append(0)
        
        alpha_int[i].append(alpha[0])

t_backstepping = back_stepping[0,:][0:100]
x_backstepping = back_stepping[1,:][0:100]


reactK_factor = Lrefstar**2*eps**2*N2star**3/Drefstar  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!


prop_c = []

print(np.array(t_induct) * Lrefstar ** 2 / Drefstar,'t_induct')

for i in range(len(k_values)):
    prop_c.append(t_induct[i]/((k_values[i]*reactK_factor)**(-0.25) / eps ))

print(prop_c[0])
prop_c = 6
t_scalings = [(k*reactK_factor)**(-0.25) / eps * prop_c for k in k_values]
# t_scalings = t_induct
# print(t_scalings)


# plt.plot(t_values-t_scalings[0],x1[0],label='k = 1e13')
# plt.plot(t_values-t_scalings[1] ,x1[1],label='k = 1e14')
# plt.plot(t_values-t_scalings[2],x1[2],label='k = 1e15')
# plt.plot(t_values-t_scalings[3],x1[3],label='k = 1e16')



plt.plot(Trefstar * t_backstepping / 3600,Lrefstar * (x_backstepping) / 1e-7,label='Asymptotic solution',linestyle='dashed')

# log-log-plot
# plt.plot(np.log(Trefstar * t_backstepping / 3600),np.log(Lrefstar * (x_backstepping) / 1e-7),label='Asymptotic solution',linestyle='dashed')


# plt.plot(t_values * Trefstar / 3600 ,(np.sqrt(2/3 * eps * D1 * (1-cs) * t_values) - D1/D3 * L3) * Lrefstar / 1e-7,label='Asymptotic solution',linestyle='dashed')
print(D1/D3*L3)
# plt.plot(np.log(t_backstepping),np.log(x_backstepping),label='Asymptotic solution')
# plt.plot(np.log(t_values),np.log(np.sqrt(2/3 * eps * D1 * (1-cs) * t_values) - D1/D3 * L3),label='Asymptotic solution')

# plt.plot(x2_nodes[0:200],alpha[0:200])
# plt.plot(x2_nodes[0:200],c2[0:200])
print(t_values[-1] * Lrefstar**2/Drefstar)
k_value_string = [r'$10^{13}$',r'$10^{14}$',r'$10^{15}$',r'$10^{16}$']
for i in range(len(k_values)):

    plt.plot((t_values-t_scalings[i]) * Trefstar / 3600,Lrefstar * np.array(x1[i]) / 1e-7,label=f'$k^*$ = {k_value_string[i]}')


    # log-log-plot
    # plt.plot(np.log((t_values-t_scalings[i]) * Trefstar / 3600),np.log(Lrefstar * np.array(x1[i]) / 1e-7),label=f'$k^*$ = {k_value_string[i]}')


    # plt.plot((t_values-t_scalings[i]) * Trefstar / 3600,Lrefstar * np.array(x1[i]) / 1e-7,label=f'$k^*$ = {k_values[i]:.2e}')

    # plt.plot(np.log(t_values-t_scalings[i]),np.log(x1[i]),label=f'$k^*$ = {k_values[i]:2e}')
    


    # m, b = np.polyfit(np.log((t_values-t_scalings[i])[90:100]),np.log(np.array(x1[i])[90:100]), 1)
    # m, b = np.polyfit(np.log((t_backstepping)[90:100]),np.log(np.array(x_backstepping)[90:100]), 1)
    # m, b = np.polyfit(np.log((t_values)[90:100]),np.log((np.sqrt(2/3 * eps * D1 * (1-cs) * t_values))[90:100]), 1)
    # print(i,m)
    
    # plt.plot(t_values * Trefstar / 3600 ,alpha_int[i],label=f'$k^*$ = {k_values[i]:.2e}')

# log-log
# plt.xlabel(r'$\mathrm{log}(t^*$/ h)')
# plt.ylabel(r'$\mathrm{log}$(Hydride thickness / $\mathrm{nm})$')

plt.xlabel(r'$t^*$/ h')
plt.ylabel(r'Hydride thickness / $\mathrm{nm}$')

plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.legend(loc='lower right')
plt.show()