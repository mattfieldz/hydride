import numpy as np
from matplotlib import pyplot as plt


L3star = 20*1.0e-7  # oxide lengthscale 10 nm = 1.e-8 m = 1.e-6 cm
L2star = 1000*1.0e-7  # bulk lengthscale  50 um = 5.e-5 m = 5.e-3 cm
Lmstar = 5000*1.0e-7



D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 @ room temp
D2star = 1.49e-10  # cm^2/s, diffusion of H in U @ room temp
D1star = 1.0e-13


D1 = D1star/D2star
D2 = 1
D3 = D3star/D2star

L3 = L3star/L2star


N1star = 4.54e-2 
N2star = 8.01e-2
N3star = 4.06e-2

k2star = 1.0e13

Castar = 1.0e-4
Csstar = 1.0e-5

Lref = L2star
Dref = D2star


epsilon = Castar / N2star


c_sol = 0.1

# non-dimensional domains
L3 = L3star / L2star
Lm = Lmstar / L2star
L2 = 1


A = epsilon / 3
t = 0
dt = -0.01
fig, axes = plt.subplots(1,2)
x1 = 0.2
x1_list = [x1]
t_list = [0]
counter = 0
while x1 > 0:
    c_i = (D1/D3 * c_sol * L3 + x1)/(x1+D1/D3 * L3)


    dxdt = - A * D1 * (c_sol-c_i)/x1

    x1 += dt * dxdt
    t += -dt
    counter += 1 
    if int(counter) == 1000:
        x1_list.append(x1)
        t_list.append(t)    
        print(x1)
        counter = 0
x1_list.reverse()
 
t_array = np.array(t_list)
x_array = np.array(x1_list)
np.savetxt(
            f"1Dexamples/backtimestepsk7.dat",
            np.array((t_array,x_array)))

print('t : ',t)
# axes[0] = plt.plot(np.log10(t_array),np.log10(x_array))

axes[0] = plt.plot((t_array)**0.5,(x_array))
plt.show()   