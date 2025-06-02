from matplotlib import pyplot as plt
import numpy as np

alpha_varying = np.loadtxt(
    f"formal_work/data1D/c_xi_0_list_alpha_varying.dat",
    )

m = 1
alpha_c = 0.98

alpha_values = alpha_varying[:,0]
c_xi_alpha = alpha_varying[:,1]

m_varying = np.loadtxt(
    f"formal_work/data1D/c_xi_0_list_m_varying.dat",
    )

cs_varying = np.loadtxt(
    f"formal_work/data1D/c_xi_0_list_cs_computed_varying.dat",
    )
cs_values = cs_varying[:,0]
c_xi_cs = cs_varying[:,1]

m_values = m_varying[:,0]
print(m_values)
c_xi_m = m_varying[:,1]

fig, axes = plt.subplots(1,1)

# axes[0].plot(m_values,c_xi_m,label=r'$\alpha_c=0.98$, $c_s=0.1$')
# axes[1].plot(alpha_values,c_xi_alpha,label=r'$m=3$, $c_s=0.1$')
# # axes[2].plot(cs_values,c_xi_cs,label=r'$\alpha_c=0.98$, $m=3$')

# axes[0].set_xlabel('m')
# axes[0].set_ylabel(r'$C_{\xi}(0)$')

# axes[1].set_xlabel(r'$\alpha_c$')
# axes[1].set_ylabel(r'$C_{\xi}(0)$')

# # axes[2].set_xlabel(r'$c_s$')
# # axes[2].set_ylabel(r'$C_{\xi}(0)$')

# axes[0].legend()
# axes[1].legend()
# # axes[2].legend()
# fig.tight_layout() 


axes.plot(cs_values,c_xi_cs,label=r'$\alpha_c=0.98$, $m=1$')

import scipy as sp

lin = sp.stats.linregress(cs_values,c_xi_cs)
print(lin)


axes.set_xlabel(r'$c_s$')
axes.set_ylabel(r'$C_{\xi}(0)$')
axes.legend()
plt.show()
