import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import scipy as sp

plt.style.use('formal_work/minorticks.mplstyle')


data = pd.read_csv('formal_work/condon_kirkpatrick/data.txt')
cs_varying = np.loadtxt(
    f"formal_work/data1D/c_xi_0_list_cs_computed_varying.dat",
    )

alpha_varying = np.loadtxt(
    f"formal_work/data1D/c_xi_0_list_alpha_new_varying.dat",
    )

print(cs_varying)

regress = sp.stats.linregress(cs_varying)


m = regress[0]
b = regress[1]
# print(regress)


arr = np.array(data)

alpha_c = 0.98
N2star = 8.01 * 1e-2

T_values = 1000 / arr[:,0]

# print('t-values',T_values)

P = 101325

R = 8.31
D_values = 1.9 * 1e-2 * np.exp(-48387/(R*T_values)) # cm^2/s
N_values = np.exp(-2.362-2305/T_values)
# P_eq = np.exp(69.32 - 14640 / T_values - 5.65 * np.log(T_values)) * np.exp(680.8/T_values) # Pa

# C_values = N2star * N_values * 4.13 * 1e-6 * np.exp(-894/T_values) * P**0.5 / (N_values + 4.13 * 1e-6 * np.exp(-894/T_values) * P**0.5) # mol/cm^3
# cs_values = N2star * N_values * 4.13 * 1e-6 * np.exp(-894/T_values) * P_eq**0.5 / (N_values * 4.13 * 1e-6 * np.exp(-894/T_values) * P_eq**0.5 ) # mol/cm^3

C_values = N2star * 4.13 * 1e-6 * np.exp(-894/T_values) * P**0.5
cs_values = 0.25 * np.exp(-52000/(R*T_values))


k_values = 1/N2star * 10.4 * np.exp(1542/T_values) # (molcm^-3)^-1 s^-1
c_xi_0_values = m * cs_values/C_values + b
i=1
# print(T_values[i],D_values[i],C_values[i],cs_values[i]/C_values[i],k_values[i])


# print(4.13 * 1e-6 * np.exp(-894/T_values) * P_eq**0.5)
# print(c_xi_0_values)
# print(m)

# 0.98/0.97 are alpha_c values

print('lol',c_xi_0_values)

V_min99 = -np.sqrt(D_values * k_values * (C_values-0*cs_values)**2 / N2star) * (-1/(3 * (1 - 0.99))) * 1e7
V_min98 = -np.sqrt(D_values * k_values * (C_values-0*cs_values)**2 / N2star) * (-1/(3 * (1 - 0.98))) * 1e7
V_min97 = -np.sqrt(D_values * k_values * (C_values-0*cs_values)**2 / N2star) * (-1/(3 * (1 - 0.97))) * 1e7
# V_max = -np.sqrt(D_values * k_values * (C_values-0*cs_values)**2 / N2star) * (-1/(3 * (1 - 0.97))) * 1e7

# print(V)

# print('cs_values',cs_values,T_values)

# print(len(V))
# print(V)
# regress_V = sp.stats.linregress(arr[0:25,0],np.log10(V[0:25]))

# print(regress_V)


# plt.scatter(arr[:,0],10**arr[:,1],marker='x',s=20,c='red')
# plt.plot(arr[:,0],(V_min99),linestyle='dashed',label=r'$\alpha_c=0.99$')
# plt.plot(arr[:,0],(V_min98),linestyle='dashed',label=r'$\alpha_c=0.98$')
# plt.plot(arr[:,0],(V_min97),linestyle='dashed',label=r'$\alpha_c=0.97$')

plt.scatter(T_values,10**arr[:,1],marker='x',s=20,c='red')
plt.plot(T_values,(V_min99),linestyle='dashed',label=r'$\alpha_c=0.99$')
plt.plot(T_values,(V_min98),linestyle='dashed',label=r'$\alpha_c=0.98$')
plt.plot(T_values,(V_min97),linestyle='dashed',label=r'$\alpha_c=0.97$')

# plt.plot(arr[:,0],np.log10(V_max))
# plt.xlabel(r'$T^* \, /\, \mathrm{K} $')
# plt.ylabel(r'$V^* / (nm\,s^{-1})$')
plt.xscale('function',functions=(lambda x : 1000/x, lambda x : 1000/x))
plt.yscale('log')


# plt.legend()
plt.savefig('formal_work/1D_code/regime2a_plots/figures_latex/regime3-data-plot.eps',format='eps',bbox_inches='tight')

plt.show()