import numpy as np
from matplotlib import pyplot as plt

plt.style.use('formal_work/minorticks.mplstyle')

T = 50 + 273.15

P = np.linspace(1e5,1e7,100)


N = np.exp(-2.362-2305/T)

base = 4.13 * 10**(-6) * np.exp(-894/T) * P**0.5

c0 = N * base / (N+base)

# plt.plot(np.log10(P),np.log10(c0))
# plt.show()

T = np.linspace(25+273.15,1000+273.15,100)

uranium_c0 = 3.32 * 10**(-7) * np.exp(-894/T) * np.sqrt(100000)
dioxide_c0 = 7.1*10**(-6) * np.exp(-12030/T) * np.sqrt(100000)

print(5.5 * 1e4 * 4.06 * 10**(-2) * (1e-6) * (10**(-5/2)))
print(100000/8.31)

# plt.plot(T,uranium_c0,label='uranium')
# plt.plot(T,dioxide_c0,label='dioxide')

plt.plot(T-273.15,np.log10(uranium_c0),label='Uranium metal')
plt.plot(T-273.15,np.log10(dioxide_c0),label='Uranium dioxide',linestyle='dashed')

plt.xlabel(r'${T\,/\,\mathrm{C^{\circ}}}$')
plt.ylabel(r'$\mathrm{log_{10}}(c_0^* \,/\, \mathrm{mol\,cm^{-3}}$)')

plt.legend()
plt.show()
