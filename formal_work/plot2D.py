from matplotlib import pyplot as plt
import numpy as np

# t_values = np.loadtxt('1Dexamples/2Ddata/tspace.dat')
t_values = [250*i+0.00001 for i in range(1,int(25))]
t_values = np.array(t_values)

n3 = 51
n2 = 401
m = 401

L3star = 10*1.0e-7  # oxide lengthscale 10 nm = 1.e-8 m = 1.e-6 cm
L2star = 1000*1.0e-7  # bulk lengthscale  1 um = 5.e-5 m = 5.e-3 cm
Lmstar = 750*1.0e-7



D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 @ room temp
D2star = 1.49e-10  # cm^2/s, diffusion of H in U @ room temp
D1star = 1.0e-13


N1star = 4.54e-2 
N2star = 8.01e-2
N3star = 4.06e-2

k2star = 2.0e13

Castar = 1.0e-4
Csstar = 1.0e-5

Lref = 1.0e-7
Dref = D2star

epsilon = Castar / N2star


c_sol = 0.1

# non-dimensional domains
L3 = L3star / Lref
Lm = Lmstar / Lref
L2 = L2star / Lref


# fig,axis = plt.subplots(2,2,subplot_kw={"projection": "3d"})
fig,axis = plt.subplots(2,2)
fig.subplots_adjust(wspace=0.5, hspace=0.5)

x2_nodes = np.loadtxt(f'formal_work/data2D/k0/no_oxide/glascott/x2_nodes_only.dat')
x3_nodes = np.loadtxt(f'formal_work/data2D/k0/1nm/x3_nodes.dat')
m_axis = np.linspace(0,Lm,m)

print(len(x2_nodes))
M, X2 = np.meshgrid(m_axis,Lref*x2_nodes)

length_vector = np.zeros(n2)

length_vector[0:n2] = np.linspace(0,L2,n2)

c2_int = []
c2_deep = []
t_c = 0
for t in t_values:
    t_c += 1
    c2 = np.loadtxt(f'formal_work/data2D/c2_{t:.2f}.dat')
    alpha = np.loadtxt(f'formal_work/data2D/alpha_{t:.2f}.dat')
    # c2_1D = np.loadtxt(f'formal_work/data1D/k1e6/c2{t:.2f}.dat')
    # alpha_1D = np.loadtxt(f'formal_work/data1D/k1e6/alpha{t:.2f}.dat')
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

    # pcm_conc = axis[0,0].plot(x2_nodes,Castar*c2[:,17],color='blue')
    # pcm_conc_1D = axis[0,0].plot(x2_nodes,Castar*c2_1D[:,1],color='red')
    
    
    # pcm_alpha = axis[0,1].plot(x2_nodes,N2star * (1-alpha[:,19]),color='blue')
    # pcm_alpha_1D = axis[0,1].plot(x2_nodes,N2star * (1-alpha_1D[:,1]),color='red')

    c2_int.append(c2[0,19])
    c2_deep.append(c2[10,19])
    
    cross_section = axis[0,1].plot(t_values[0:t_c],c2_int)
    cross_section_deeper = axis[0,1].plot(t_values[0:t_c],c2_deep)

    print(np.shape(c2))
    wireframe_conc = axis[1,0].contour(M,X2,c2)
    # wirefram_alpha = axis[1,1].contour(M,X2,1-alpha)
    print(c2[:,0])
    # pcm_beta = axis[1,0].imshow(c2*Castar,aspect='auto',extent=[0,Lmstar * 1e7,-L2star * 1e7,0])
    # pcm_gamma = axis[1,1].imshow(c2[]*N2star,aspect='auto',extent=[0,Lmstar * 1e7,-L2star * 1e7,0])
    fig.suptitle(('dim t',t*Lref**2/(Dref),'non-dim: ',t))
    
    
    plt.pause(1)
    