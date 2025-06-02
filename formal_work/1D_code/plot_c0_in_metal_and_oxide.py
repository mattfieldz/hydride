from matplotlib import pyplot as plt
import numpy as np

R = 8.31 * 10**(-3)

P = 100000

T = np.linspace(273.15,273.15+200,100)

c_o = 5.5 * 1e4 * np.exp(-100/(R*T)) * np.sqrt(P/(1e6)) * 1e6 * 4*1e-2
c_m = 8.03 * 1e-2 * 4.13 * 1e-6 * np.exp(-894/T) * np.sqrt(P) 

D_o = 0.037 * np.exp(-60/(R*T))


plt.plot(np.log10(T-273.15),np.log10(c_o))
plt.plot(np.log10(T-273.15),np.log10(c_m))
# plt.plot(T-273.15,D_o)
plt.show()