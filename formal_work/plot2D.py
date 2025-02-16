from matplotlib import pyplot as plt
import numpy as np

# t_values = np.loadtxt('1Dexamples/2Ddata/tspace.dat')
t_values = [2*i+0.00001 for i in range(1,int(440))]
t_values = np.array(t_values)

n3 = 51
n2 = 101
m = 40

L3star = 10*1.0e-7  # oxide lengthscale 10 nm = 1.e-8 m = 1.e-6 cm
L2star = 1000*1.0e-7  # bulk lengthscale  1 um = 5.e-5 m = 5.e-3 cm
Lmstar = 1e3*1.0e-7



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


fig,axis = plt.subplots(2,2,subplot_kw={"projection": "3d"})
# fig,axis = plt.subplots(2,2)
fig.subplots_adjust(wspace=0.5, hspace=0.5)

x2_nodes = np.loadtxt(f'formal_work/data2D/k1e6/De-15/x2_nodes.dat')
x3_nodes = np.loadtxt(f'formal_work/data2D/k1e6/De-15/x3_nodes.dat')
m_axis = np.linspace(0,Lm,m)


X2, M = np.meshgrid(x2_nodes,m_axis)

length_vector = np.zeros(n2)

length_vector[0:n2] = np.linspace(0,L2,n2)



for t in t_values:
    c2 = np.loadtxt(f'formal_work/data2D/k1e6/De-15/c2{t:.2f}.dat')
    alpha = np.loadtxt(f'formal_work/data2D/k1e6/De-15/alpha{t:.2f}.dat')
    c2_1D = np.loadtxt(f'formal_work/data1D/k1e6/c2{t:.2f}.dat')
    alpha_1D = np.loadtxt(f'formal_work/data1D/k1e6/alpha{t:.2f}.dat')
    axis[0,0].cla()
    axis[0,1].cla()
    axis[1,0].cla()
    axis[1,1].cla()

    axis[0,0].set_title('[H] concentration cross section')
    axis[0,1].set_title('[U] concentration cross section')
    axis[1,0].set_title('[H] / $\mathrm{mol}$ $\mathrm{cm}{-3}$')
    axis[1,1].set_title('[UH3] / $\mathrm{mol}$ $\mathrm{cm}^{-3}$')

    axis[1,0].set_xlabel('$x / \mathrm{nm}$')
    axis[1,1].set_xlabel('$x / \mathrm{nm}$')

    axis[1,0].set_ylabel('$y / \mathrm{nm}$')
    axis[1,1].set_ylabel('$y / \mathrm{nm}$')

    axis[0,0].set_xlabel('$y / \mathrm{nm}$')
    axis[0,0].set_ylabel('[H] / $\mathrm{mol}$ $\mathrm{cm}{-3}$')


    axis[0,1].set_xlabel('$y / \mathrm{nm}$')
    axis[0,1].set_ylabel('[U] / $\mathrm{mol}$ $\mathrm{cm}{-3}$')

    pcm_conc = axis[0,0].plot(x2_nodes,Castar*c2[:,0],color='blue')
    pcm_conc_1D = axis[0,0].plot(x2_nodes,Castar*c2_1D[:,1],color='red')
    
    
    pcm_alpha = axis[0,1].plot(x2_nodes,N2star * (1-alpha[:,0]),color='blue')
    pcm_alpha_1D = axis[0,1].plot(x2_nodes,N2star * (1-alpha_1D[:,1]),color='red')
    
    wireframe_conc = axis[1,0].plot_surface(X2,M,c2.T)
    wirefram_alpha = axis[1,1].plot_surface(X2,M,1-alpha.T)
    print(c2[:,0])
    # pcm_beta = axis[1,0].imshow(c2[0::2]*Castar,aspect='auto',extent=[0,Lmstar * 1e7,-L2star * 1e7,0])
    # pcm_gamma = axis[1,1].imshow(c2[1::2]*N2star,aspect='auto',extent=[0,Lmstar * 1e7,-L2star * 1e7,0])
    fig.suptitle(('dim t',t*Lref**2/(Dref),'non-dim: ',t))
    
    plt.pause(0.01)