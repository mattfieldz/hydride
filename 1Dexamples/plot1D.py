from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
t_values = np.loadtxt('1Dexamples/1Ddata/tspace.dat')


e1 = 2
n3 = 51
n2 = 101


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
length_vector = np.zeros(n2+n3)
length_vector[0:n3] = np.linspace(-L3,0,n3)
length_vector[n3:n2+n3] = np.linspace(0,L2,n2)

x = np.zeros(e1*(n2+n3))
x1 = []
x2 = []
A = D1star * epsilon / 3


t_values_dim = t_values * Lref**2 / Dref



c2_interface = []
t_c = 0
for t in t_values:
    t_c += 1
    c2 = np.loadtxt(f'1Dexamples/1Ddata/c2_t{t}')
    c3 = np.loadtxt(f'1Dexamples/1Ddata/c3_t{t}')

    x[0:e1*n3] = c3
    x[e1*n3:e1*(n2+n3)] = c2

    c2_interface.append(c2[0])



    spline = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0,L2,n2),(c2[0::2]-c_sol)**3 * (1-c2[1::2]),k=4)
    

    cr_points = spline.derivative().roots()
    
    cr_values = spline(np.append(cr_points,(0,L2)))
    cr_maximal_index = np.argmax(cr_values)

    if len(cr_points) > 0:
        x1.append(np.append(cr_points,(0,L2))[cr_maximal_index])
    else:
        x1.append(0)
       

    c2h = np.loadtxt(f'1Dexamples/1Ddata/c2_th{t}')
    c3h = np.loadtxt(f'1Dexamples/1Ddata/c3_th{t}')

    x[0:e1*n3] = c3h
    x[e1*n3:e1*(n2+n3)] = c2h

    # c2_interface.append(c2h[0])



    splineh = sp.interpolate.InterpolatedUnivariateSpline(np.linspace(0,L2,n2),(c2h[0::2]-c_sol)**3 * (1-c2h[1::2]),k=4)
    

    cr_pointsh = splineh.derivative().roots()
    
    cr_valuesh = splineh(np.append(cr_pointsh,(0,L2)))
    cr_maximal_indexh = np.argmax(cr_valuesh)

    if len(cr_pointsh) > 0:
        x2.append(np.append(cr_pointsh,(0,L2))[cr_maximal_indexh])
    else:
        x2.append(0)
  
    if t_c > 20:
        m, b = np.polyfit(np.log10(t_values_dim[t_c-5:t_c]),np.log10(np.array(x1)[t_c-5:t_c]*L2star-D1star/D3star*L3star), 1)
        print(m,b)
        # print(np.log10(t_values_dim[70:t_c],),np.log10(np.array(x1)[70:t_c]*L2star))
    else:
        m,b=0,0

    axis[0,0].cla()
    axis[0,1].cla()
    axis[1,0].cla()
    axis[1,1].cla()

    axis[0,0].set_title('H concentration cross section')
    axis[0,1].set_title('UH3 concentration cross section')
    axis[1,0].set_title('H conc ')
    axis[1,1].set_title('UH3 conc')


    pcm_conc = axis[0,0].plot(length_vector,x[0::2])
    
    
    pcm_alpha = axis[0,1].plot(length_vector,1-x[1::2])

    pcm_beta = axis[1,0].plot(t_values[0:t_c],c2_interface)

    # pcm_gamma = axis[1,1].plot(c2[0::2]**3 * (1-c2[1::2]))
    # pcm_gamma = axis[1,1].plot(t_values_dim[0:t_c],np.array(x1)*L2star)
    plot_A = axis[1,1].plot(t_values_dim[0:t_c],np.sqrt(2*A*(1-c_sol)*t_values_dim[0:t_c]))
    # plot_log = axis[1,1].plot(np.log10(t_values_dim[0:t_c],),np.log10(np.array(x1)*L2star))
    # plot_log = axis[1,1].plot(np.log10(t_values_dim[0:t_c],),np.log10(np.array(x2)*L2star*np.sqrt(2)))

    plot_thalf = axis[1,1].plot(t_values[0:t_c],np.array(x1),label='epsilon 0.001')
    plot_thalf2 = axis[1,1].plot(t_values[0:t_c],np.array(x2)*np.sqrt(2),label='epsilon 0.0005')
    axis[1,1].legend()
    # pcm_gamma = axis[1,1].plot(c2[1::2])

    fig.suptitle((t*Lref**2/(Dref),'non-dim: ',t))

    plt.pause(0.01)
    