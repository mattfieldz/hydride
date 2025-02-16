from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
t_values = [1*i+0.00001 for i in range(1,int(440))]
t_values = np.array(t_values)



def x_from_X(X):
    """Defines the (non-uniformly distributed) physical coordinate 'x'
    in terms of the (uniformly distributed) computational coordinate
    'X'. Returns: not only x(X) but also x'(X) and x''(X).

    """
    BX = 0.01  # smaller BX leads to more NON-uniformity
    # +0*X below prevents numpy from returning a scalar when acting
    # on an iterable container for X
    return (X + BX) * (X + BX) - BX * BX, 2 * (X + BX)



e1 = 2
n3 = 101
n2 = 2001


L3star = 10*1.0e-7  # oxide lengthscale 10 nm = 1.e-8 m = 1.e-6 cm

L2star = 1000*1.0e-7  # bulk lengthorscale  50 um = 5.e-5 m = 5.e-3 cm
Lmstar = 5000*1.0e-7



D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 @ room temp
D2star = 1.49e-10  # cm^2/s, diffusion of H in U @ room temp
D1star = 1.0e-12


D1 = D1star/D2star
D3 = D3star/D2star

L3 = L3star/L2star


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

x = np.zeros(e1*(n2+n3))
x1 = []
x2 = []
A = D1star/Dref * epsilon / 3


t_values_dim = t_values * Lref**2 / Dref

k=1e7

c2_interface = []
alpha_interface = []
left_side_list = []
right_side_list = []
c_i_asymptotic_list = []
c_i_list = []
t_c = 0
D = []
eps = np.exp(-epsilon*(D3/D1)*t_values*(1-c_sol)**3)
print(eps)

backtimesteps =  np.loadtxt(
            f"1Dexamples/backtimestepsk7.dat")

spline_back = sp.interpolate.InterpolatedUnivariateSpline(backtimesteps[0],backtimesteps[1],k=4)
    
back_steps = []
back_steps_time = []

for t in t_values:
    t_c += 1
    
    c3_array =  np.loadtxt(
            f"1Dexamples/datarich/c3_D112_k1e7{t:.2f}.dat")
    c2_array =  np.loadtxt(
            f"1Dexamples/datarich/c2_D112_k1e7{t:.2f}.dat")
    alpha_array =  np.loadtxt(
            f"1Dexamples/datarich/alpha_D112_k1e7{t:.2f}.dat")
    
    # c3_array =  np.loadtxt(
    #         f"1Dexamples/datarich/c3_k1e5Deq{t:.2f}.dat")
    # c2_array =  np.loadtxt(
    #         f"1Dexamples/datarich/c2_k1e5Deq{t:.2f}.dat")
    # alpha_array =  np.loadtxt(
    #         f"1Dexamples/datarich/alpha_k1e5Deq{t:.2f}.dat")

    x3_nodes = c3_array[:,0]
    c3 = c3_array[:,1]
    x2_nodes = c2_array[:,0]
    c2 = c2_array[:,1]
    alpha = alpha_array[:,1]
    
    c2_interface.append(c2[0])
    alpha_interface.append(alpha[0])

    spline = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,(c2-c_sol)**3*(alpha),k=4)
    

    cr_points = spline.derivative().roots()
    
    cr_values = spline(np.append(cr_points,(0,L2)))
    cr_maximal_index = np.argmax(cr_values)

    # if len(cr_points) > 0:
    #     x1.append(np.append(cr_points,(0,L2))[cr_maximal_index])
    # else:
    #     x1.append(0)

    spline = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,alpha-0.0001)
    cr_points = spline.roots()
    if len(cr_points) > 0:
        x1.append(cr_points[0])
        # print(x1,'lol')
    else:
        x1.append(0)

    if t_c > 100:
        m, b = np.polyfit(np.log10(t_values[t_c-50:t_c]),np.log10(np.array(x1)[t_c-50:t_c]), 1)
        print(m,b)
        # print(np.log10(t_values_dim[70:t_c],),np.log10(np.array(x1)[70:t_c]*L2star))
    else:
        m,b=0,0


    c_i = c2[0]
    c_i_list.append(c_i)
    left_side = D3 * (c_i-1)/L3
    right_side = D1 * (c_sol - c_i)/x1[-1]

    c_i_asymptotic = (D1/D3 * c_sol * L3 + x1[-1])/(x1[-1] + D1/D3*L3)
    # print(c_i,c_i_asymptotic)
    # print(c2[100]-c_sol)
    left_side_list.append(left_side)

    right_side_list.append(right_side)    

    c_i_asymptotic_list.append(c_i_asymptotic)

    back_steps.append(spline_back(t))


    # axis[0,0].cla()
    # axis[0,1].cla()
    # axis[1,0].cla()
    # axis[1,1].cla()

    # axis[0,0].set_title('H concentration at interface',fontsize=10,y=1)
    # axis[0,1].set_title('U concentration',fontsize=10)
    # axis[1,0].set_title('slope comparison',fontsize=10)
    # axis[1,1].set_title('Hydride thickness',fontsize=10)

    # axis[0,0].set_xlabel('')
    # axis[0,0].set_ylabel('c')

    # axis[0,1].set_xlabel('')
    # axis[0,1].set_ylabel('[U]')


    # # axis[1,0].set_xlabel('$t^*/s$')
    # # axis[1,0].set_ylabel('[H] / $\mathrm{molcm}^{-3}$')

    # axis[1,0].set_xlabel('t')
    # axis[1,0].set_ylabel('gradient approximation')

    # axis[1,1].set_xlabel('$t$')
    # axis[1,1].set_ylabel('x1')


    # # pcm_conc = axis[0,0].plot(np.append(x3_nodes,x2_nodes),np.append(c3,c2))
    
    
    # pcm_alpha = axis[0,1].plot(x2_nodes,alpha)
    # plot_single_point = axis[0,1].plot(x1,np.ones(len(x1)))
    # plot_single_point2 = axis[0,1].plot(np.sqrt(2*(D1)*A*(1-c_sol)*t_values[0:t_c]),0.5*np.ones(len(x1)))

    # # pcm_beta = axis[1,0].plot(t_values_dim[0:t_c],c2_interface)



    # pcm_gamma = axis[0,0].plot(t_values[0:t_c],c2_interface,label='numerics')
    # plotting_asymptotic_c_interface = axis[0,0].plot(t_values[0:t_c],c_i_asymptotic_list,label='asymptotic')

    # axis[0,0].legend()
    # # plot_react = axis[1,0].plot(x2_nodes,(c2-c_sol))
    
    # plot_sides = axis[1,0].plot(t_values[0:t_c],left_side_list,label='slope in oxide')
    # plot_right_side = axis[1,0].plot(t_values[0:t_c],right_side_list,label='slope in hydride')

    # axis[1,0].legend()

    # # pcm_gamma = axis[1,1].plot(t_values_dim[0:t_c],np.array(x1)*L2star)
    # # plot_A = axis[1,1].plot(np.log10(t_values_dim[0:t_c]),np.log10(np.sqrt(2*A*(1-c_sol)*t_values_dim[0:t_c])-epsilon*L3star))
    # # plot_log = axis[1,1].plot(np.log10(t_values_dim[0:t_c],),np.log10(np.array(x1)*L2star))
    # # plot_line_log = axis[1,1].plot(np.log10(t_values_dim[0:t_c]),m*np.log10(t_values_dim[0:t_c])+b,label=f'm={m:.2f},b={b:.2f}')    

    # # plot_alpha = axis[1,1].plot(t_values_dim[0:t_c], alpha_interface)
    
    # spline_alpha = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,alpha,k=4)
    # spline_conc = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,c2,k=4)
    
    # # print(x1[-1])
    # # print(spline_alpha(x1[-1]))
    # # print(spline_conc(x1[-1]))
    
    
    # # plot_log = axis[1,1].plot(np.log10(t_values_dim[0:t_c],),np.log10(np.array(x2)*L2star*np.sqrt(2)))
    
    # plot_thalf = axis[1,1].plot(t_values[0:t_c],np.array(x1),label='Numerics')

    
    # plot_back = axis[1,1].plot(t_values[0:t_c],back_steps,label='backsteps')
    # plot_A = axis[1,1].plot((t_values[0:t_c]),(np.sqrt(2*(D1)*epsilon/3*(1-c_sol)*((t_values[0:t_c])))),label='Approx')
    # # plot_A = axis[1,1].plot((t_values[0:t_c]+200),(2 * np.array(right_side_list) * epsilon/3 * (t_values[0:t_c])**0.5),label='Approx')
    
    
    # # plot_A = axis[1,1].plot((t_values[0:t_c]),(np.sqrt(-(2*epsilon*(c_sol-1)*(D1 - (D1-1)*eps[0:t_c]))/(3*(1-eps[0:t_c]))*((t_values[0:t_c])))),label='Approx 2')
    
    # # plot_thalf2 = axis[1,1].plot(t_values_dim[0:t_c],np.array(x2)*L2star*np.sqrt(2),label='epsilon 0.0005')
    # axis[1,1].legend()

    # # pcm_gamma = axis[1,1].plot(c2[1::2])

    # fig.suptitle((t*Lref**2/(Dref),'non-dim: ',t))

    # plt.pause(0.01)




axis[0,0].set_title('H concentration at interface',fontsize=10,y=1)
axis[0,1].set_title('U concentration',fontsize=10)
axis[1,0].set_title('slope comparison',fontsize=10)
axis[1,1].set_title('Hydride thickness',fontsize=10)

axis[0,0].set_xlabel('')
axis[0,0].set_ylabel('c')

axis[0,1].set_xlabel('')
axis[0,1].set_ylabel('[U]')


# axis[1,0].set_xlabel('$t^*/s$')
# axis[1,0].set_ylabel('[H] / $\mathrm{molcm}^{-3}$')

axis[1,0].set_xlabel('t')
axis[1,0].set_ylabel('gradient approximation')

axis[1,1].set_xlabel('$t$')
axis[1,1].set_ylabel('x1')


# pcm_conc = axis[0,0].plot(np.append(x3_nodes,x2_nodes),np.append(c3,c2))


pcm_alpha = axis[0,1].plot(x2_nodes,alpha)
plot_single_point = axis[0,1].plot(x1,np.ones(len(x1)))
plot_single_point2 = axis[0,1].plot(np.sqrt(2*(D1)*A*(1-c_sol)*t_values[0:t_c]),0.5*np.ones(len(x1)))

# pcm_beta = axis[1,0].plot(t_values_dim[0:t_c],c2_interface)



pcm_gamma = axis[0,0].plot(t_values[0:t_c],c2_interface,label='numerics')
plotting_asymptotic_c_interface = axis[0,0].plot(t_values[0:t_c],c_i_asymptotic_list,label='asymptotic')

axis[0,0].legend()
# plot_react = axis[1,0].plot(x2_nodes,(c2-c_sol))

plot_sides = axis[1,0].plot(t_values[0:t_c],left_side_list,label='slope in oxide')
plot_right_side = axis[1,0].plot(t_values[0:t_c],right_side_list,label='slope in hydride')

axis[1,0].legend()

# pcm_gamma = axis[1,1].plot(t_values_dim[0:t_c],np.array(x1)*L2star)
# plot_A = axis[1,1].plot(np.log10(t_values_dim[0:t_c]),np.log10(np.sqrt(2*A*(1-c_sol)*t_values_dim[0:t_c])-epsilon*L3star))
# plot_log = axis[1,1].plot(np.log10(t_values_dim[0:t_c],),np.log10(np.array(x1)*L2star))
# plot_line_log = axis[1,1].plot(np.log10(t_values_dim[0:t_c]),m*np.log10(t_values_dim[0:t_c])+b,label=f'm={m:.2f},b={b:.2f}')    

# plot_alpha = axis[1,1].plot(t_values_dim[0:t_c], alpha_interface)

spline_alpha = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,alpha,k=4)
spline_conc = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,c2,k=4)

# print(x1[-1])
# print(spline_alpha(x1[-1]))
# print(spline_conc(x1[-1]))


# plot_log = axis[1,1].plot(np.log10(t_values_dim[0:t_c],),np.log10(np.array(x2)*L2star*np.sqrt(2)))

plot_thalf = axis[1,1].plot(t_values[0:t_c],np.array(x1)[0:t_c],label='Numerics')


plot_back = axis[1,1].plot(t_values[0:t_c-42]+42,np.array(back_steps)[0:t_c-42],label='backsteps')
plot_A = axis[1,1].plot((t_values[0:t_c-42]+42),(np.sqrt(2*(D1)*epsilon/3*(1-c_sol)*((t_values[0:t_c-42])))),label='Approx')
# plot_A = axis[1,1].plot((t_values[0:t_c]+200),(2 * np.array(right_side_list) * epsilon/3 * (t_values[0:t_c])**0.5),label='Approx')


# plot_A = axis[1,1].plot((t_values[0:t_c]),(np.sqrt(-(2*epsilon*(c_sol-1)*(D1 - (D1-1)*eps[0:t_c]))/(3*(1-eps[0:t_c]))*((t_values[0:t_c])))),label='Approx 2')

# plot_thalf2 = axis[1,1].plot(t_values_dim[0:t_c],np.array(x2)*L2star*np.sqrt(2),label='epsilon 0.0005')
axis[1,1].legend()

# pcm_gamma = axis[1,1].plot(c2[1::2])

fig.suptitle((t*Lref**2/(Dref),'non-dim: ',t))

plt.show()

