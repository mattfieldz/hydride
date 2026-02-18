import numpy as np
from matplotlib import pyplot as plt

# reaction order
m = 1
# temperature in celsius, pressure in atm
tempC = 25
pressure = 1

tempK = 273.15 + tempC

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 1e2 * 1e3 * 1.0e-7  # bulk domain lengthscale 1000nm=1um
L3star = 10 * 1.0e-7  # oxide domain lengthscale 10nm
D1star = 1.0e-13  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 1.94e-2 * np.exp(-5570/tempK)

D3star = 0.037 * np.exp(-7200/tempK) # cm^2/s, diffusion of H in UO2

N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 0.2446 * np.exp(-52000 / (8.31 * tempK)) # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castarck = N2star * 4.13e-6 * np.exp(-894/tempK) * np.sqrt(1e5 * pressure)  # mol/cm^3, surface value for [H] -- *ad-hoc* value
# Castar = 1e-1 * np.exp(-4000/tempK)


kstar = 10.4/N2star * np.exp(1542/tempK) # reaction constant in cm^3/mol/s

print(kstar)
print(D2star)
print(Castarck)


alpha_c = 0.98

print('velocity',1e7*np.sqrt(2*D2star*kstar*(Castarck-Csstar)**(m+1)/(3*N2star*(m+1)))/(1-alpha_c))

