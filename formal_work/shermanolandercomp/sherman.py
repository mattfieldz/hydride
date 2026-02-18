import numpy as np
from matplotlib import pyplot as plt

R = 8.31

P = np.array([0.01,0.1,1])
T = np.linspace(298,298+200,10)

N_density_UO2 = 10.97
# Molar_mass = 270.03
Molar_mass = 1

S_atomic = 5.5 * 10**4 * np.exp(-100000/(R*T)) # atomic ppm

S_mol = N_density_UO2 / (Molar_mass * 10**6) * S_atomic

print(S_mol)


