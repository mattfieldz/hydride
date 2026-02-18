import numpy as np

mu = 33.6
n_mh = 0.06098
b = 2.328
sigma_f = 2.8
gamma = 14.7

f_p = 0.64
R0 = 30
V_c = 0.15
Rp = 1 * 1e1 * 1e3  # 1 um

D = 1/ (2*gamma + 0.5 * mu * b)

A = 18.573
B = -7.100
C = 5.375

print(D)

a = C*D
bb = 2/3 * D**(2/3) * B  
c = 1/3 * D**(1/3) * A 


print(bb**3-4*a*c)


x_sol = -bb + np.sqrt(bb**2 - 4*a*c) / (2*a)

print(x_sol)

print(1/np.pi * (3*f_p/(4*np.pi))**(-2/3))