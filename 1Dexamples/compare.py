import numpy as np
from matplotlib import pyplot as plt
import scipy as sp

values = [1e2,1e3,1e4,1e5]
colors = ['red','blue','green','purple']
fig, axes = plt.subplots(1,2)

Castar = 1e-4
N2star = 8.01e-2


axes[0].set_title('Hydrogen concentration at interface')
axes[1].set_title('Uranium concentration at interface')

axes[0].set_xlabel('$t^*/s$')
axes[1].set_xlabel('$t^*/s$')

axes[0].set_ylabel('$[H]/\mathrm{mol}$ $\mathrm{cm}^{-3}$')
axes[1].set_ylabel('$[U]/\mathrm{mol}$ $\mathrm{cm}^{-3}$')

count = 0
for i in values:
    
    asym_t = np.loadtxt(f'1Dexamples/k{i}/asympttime.dat')
    asym_h = np.loadtxt(f'1Dexamples/k{i}/asymptuh3.dat')
    asym_hconc = np.loadtxt(f'1Dexamples/k{i}/asymptH.dat')


    num_t = np.loadtxt(f'1Dexamples/k{i}/fullnumtimes.dat')
    num_h = np.loadtxt(f'1Dexamples/k{i}/fullnumuh3.dat')
    num_hconc = np.loadtxt(f'1Dexamples/k{i}/fullnumH.dat')

    asym_spline = sp.interpolate.InterpolatedUnivariateSpline(asym_t,asym_hconc)
    num_spline = sp.interpolate.InterpolatedUnivariateSpline(num_t,num_hconc)

    asym_rescale = asym_spline.__call__(np.linspace(500,2000,2000))
    num_rescale = num_spline.__call__(np.linspace(500,2000,2000))
    
    

    mean_error = np.sum((asym_rescale-num_rescale)/num_rescale)/2000
    print(mean_error)

    axes[0].plot(asym_t[0:1000],Castar*asym_hconc[0:1000],linestyle='dashed',color=colors[count])
    axes[0].plot(num_t[0:5000],Castar*num_hconc[0:5000],color=colors[count])

    axes[1].plot(asym_t,N2star*asym_h,linestyle='dashed',color=colors[count])
    axes[1].plot(num_t,N2star*num_h,color=colors[count],label=f'k = {i}')
    axes[1].legend(loc="upper right")
    count += 1
    print(count)
plt.show()