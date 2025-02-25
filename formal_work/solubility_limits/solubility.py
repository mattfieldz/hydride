import numpy as np

E_a = -388
A = -0.688
P = 10
T = 298


S = np.sqrt(P) * np.exp(A+E_a/T)
print(S)

# convert to molarity
# Loui suggests that the solubility is orders of magnitude larger for reaction as need some ratio with mass of uranium