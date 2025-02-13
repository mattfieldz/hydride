from matplotlib import pyplot as plt
import numpy as np
from matplotlib import animation
t_values = np.loadtxt('1Dexamples/2Ddata/tspace.dat')


n3 = 51
n2 = 101
m = 20

L3star = 10*1.0e-7  # oxide lengthscale 10 nm = 1.e-8 m = 1.e-6 cm
L2star = 1000*1.0e-7  # bulk lengthscale  50 um = 5.e-5 m = 5.e-3 cm
Lmstar = 5000*1.0e-7



D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 @ room temp
D2star = 1.49e-10  # cm^2/s, diffusion of H in U @ room temp
D1star = 1.0e-13


N1star = 4.54e-2 
N2star = 8.01e-2
N3star = 4.06e-2

k2star = 2.0e13

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

fig,axis = plt.subplots(2,2)


c2 = np.loadtxt(f'1Dexamples/2Ddata/c2_t{t_values[0]}')
pcm_beta = axis[1,0].imshow(c2[0::2],aspect='auto',extent=[-Lm/2,Lm/2,-L2,0])
pcm_gamma = axis[1,1].imshow(c2[1::2],aspect='auto',extent=[-Lm/2,Lm/2,-L2,0])
cb1 = fig.colorbar(pcm_beta,ax=axis[1,0])



length_vector = np.zeros(n2)

length_vector[0:n2] = np.linspace(0,L2,n2)
def update(frame):
    t = t_values[frame]
    c2 = np.loadtxt(f'1Dexamples/2Ddata/c2_t{t}')
    c2_1D = np.loadtxt(f'1Dexamples/1Ddata/c2_t{t}')
    axis[0,0].cla()
    axis[0,1].cla()
    axis[1,0].cla()
    axis[1,1].cla()
    fig.cla()
    fig,axis = plt.subplots(2,2)
    


    axis[0,0].set_title('H concentration cross section')
    axis[0,1].set_title('UH3 concentration cross section')
    axis[1,0].set_title('H conc ')
    axis[1,1].set_title('UH3 conc')


    pcm_conc = axis[0,0].plot(length_vector,c2[0::2,int(m/2)],color='blue')
    pcm_conc_1D = axis[0,0].plot(length_vector,c2_1D[0::2],color='red')
    
    
    pcm_alpha = axis[0,1].plot(length_vector,1-c2[1::2,int(m/2)],color='blue')
    pcm_alpha_1D = axis[0,1].plot(length_vector,1-c2_1D[1::2],color='red')

    pcm_beta = axis[1,0].imshow(c2[0::2],aspect='auto',extent=[-Lm/2,Lm/2,-L2,0])
    pcm_gamma = axis[1,1].imshow(c2[1::2],aspect='auto',extent=[-Lm/2,Lm/2,-L2,0])
    fig.suptitle(('dim t',t*Lref**2/(Dref),'non-dim: ',t))
    cb1 = fig.colorbar(pcm_beta,ax=axis[1,0])
    cb2 = fig.colorbar(pcm_gamma,ax=axis[1,1])

    return (pcm_beta)
ani = animation.FuncAnimation(fig=fig, func=update,frames=len(t_values),interval=30)
plt.show()