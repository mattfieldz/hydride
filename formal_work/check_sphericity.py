import numpy as np
from matplotlib import pyplot as plt

x2_nodes = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/x2_nodes.dat')



c2 = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/c2_401{250:.2f}.dat')
c2 = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/c2_201{14000:.2f}.dat')

plt.plot(x2_nodes[0:150],c2[:,0][0:150])
plt.plot(np.linspace(0,x2_nodes[-1],201)[0:50],c2[0,0:50],label='radial')
plt.legend()
plt.show()