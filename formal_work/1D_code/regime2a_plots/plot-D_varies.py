
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
plt.style.use('formal_work/minorticks.mplstyle')

plt.margins(0)


# tvalues = np.loadtxt(f"formal_work/data1D/lambdaP/tvalues_vary_diff_cs0.dat")
tvalues = np.loadtxt(f"formal_work/data1D/lambdaP/tvalues_vary_diff_cs0.dat")


print(len(tvalues))



hydride_interface = []
hydride_interface_D1e0 = []
hydride_interface_D1e1 = []
hydride_interface_D1e2 = []
hydride_interface_D1e3 = []


hydride_interface_mu10 = []
hydride_interface_mu100 = []

c0_list = []
alpha0_list = []

P = 1
# threshold_alpha = 0.99 - 0.1 * P
threshold_alpha = 0.99
for t in tvalues:
    alpha_tot = np.loadtxt(
                f"formal_work/data1D/lambdaP/alpha_condon_oxide_vary_diff_cs0{t:.7f}.dat")
    c2_tot = np.loadtxt(
            f"formal_work/data1D/lambdaP/c2_condon_oxide_vary_diff_cs0{t:.7f}.dat")
    
    alpha_tot_D1e3 = np.loadtxt(
                f"formal_work/data1D/lambdaP/alpha_condon_oxide_vary_diff_cs0_Dinfty1e3{t:.7f}.dat")
    c2_tot_D1e3 = np.loadtxt(
            f"formal_work/data1D/lambdaP/c2_condon_oxide_vary_diff_cs0_Dinfty1e3{t:.7f}.dat")

    c2 = c2_tot[:,1]
    alpha = alpha_tot[:,1]
    x2_nodes = alpha_tot[:,0]
    
    
    


    interpolated_alpha = sp.interpolate.CubicSpline(x2_nodes, alpha - threshold_alpha)
    roots = interpolated_alpha.roots(extrapolate=False)
    if len(roots) > 0:
        hydride_interface.append(roots[0])
    else:
        hydride_interface.append(0)
    alpha_tot_D1e0 = np.loadtxt(
                f"formal_work/data1D/lambdaP/alpha_condon_oxide_vary_diff_cs0_Dinfty1e0{t:.7f}.dat")
    alpha_D1e0 = alpha_tot_D1e0[:,1]
    interpolated_alpha_D1e0 = sp.interpolate.CubicSpline(x2_nodes, alpha_D1e0 - threshold_alpha)
    roots_D1e0 = interpolated_alpha_D1e0.roots(extrapolate=False)
    if len(roots_D1e0) > 0:
        hydride_interface_D1e0.append(roots_D1e0[0])
    else:
        hydride_interface_D1e0.append(0)

    alpha_tot_D1e1 = np.loadtxt(
                f"formal_work/data1D/lambdaP/alpha_condon_oxide_vary_diff_cs0_Dinfty1e1{t:.7f}.dat")
    alpha_D1e1 = alpha_tot_D1e1[:,1]
    interpolated_alpha_D1e1 = sp.interpolate.CubicSpline(x2_nodes, alpha_D1e1 - threshold_alpha)
    roots_D1e1 = interpolated_alpha_D1e1.roots(extrapolate=False)
    if len(roots_D1e1) > 0:
        hydride_interface_D1e1.append(roots_D1e1[0])
    else:
        hydride_interface_D1e1.append(0)

    alpha_tot_D1e2 = np.loadtxt(
                f"formal_work/data1D/lambdaP/alpha_condon_oxide_vary_diff_cs0_Dinfty1e2{t:.7f}.dat")
    alpha_D1e2 = alpha_tot_D1e2[:,1]
    interpolated_alpha_D1e2 = sp.interpolate.CubicSpline(x2_nodes, alpha_D1e2 - threshold_alpha)
    roots_D1e2 = interpolated_alpha_D1e2.roots(extrapolate=False)
    if len(roots_D1e2) > 0:
        hydride_interface_D1e2.append(roots_D1e2[0])
    else:
        hydride_interface_D1e2.append(0)


    alpha_tot_mu10 = np.loadtxt(
                f"formal_work/data1D/lambdaP/alpha_condon_oxide_vary_diff_cs0_mu10{t:.7f}.dat")
    alpha_mu10 = alpha_tot_mu10[:,1]
    interpolated_alpha_mu10 = sp.interpolate.CubicSpline(x2_nodes, alpha_mu10 - threshold_alpha)
    roots_mu10 = interpolated_alpha_mu10.roots(extrapolate=False)
    if len(roots_mu10) > 0:
        hydride_interface_mu10.append(roots_mu10[0])
    else:
        hydride_interface_mu10.append(0)

    alpha_tot_mu100 = np.loadtxt(
                f"formal_work/data1D/lambdaP/alpha_condon_oxide_vary_diff_cs0_mu100{t:.7f}.dat")
    alpha_mu100 = alpha_tot_mu100[:,1]
    interpolated_alpha_mu100 = sp.interpolate.CubicSpline(x2_nodes, alpha_mu100 - threshold_alpha)
    roots_mu100 = interpolated_alpha_mu100.roots(extrapolate=False)
    if len(roots_mu100) > 0:
        hydride_interface_mu100.append(roots_mu100[0])
    else:
        hydride_interface_mu100.append(0)


    alpha_D1e3 = alpha_tot_D1e3[:,1]
    interpolated_alpha_D1e3 = sp.interpolate.CubicSpline(x2_nodes, alpha_D1e3 - threshold_alpha)
    roots_D1e3 = interpolated_alpha_D1e3.roots(extrapolate=False)
    if len(roots_D1e3) > 0:
        hydride_interface_D1e3.append(roots_D1e3[0])
    else:
        hydride_interface_D1e3.append(0)

    c0_list.append(c2[0])
    alpha0_list.append(alpha[0])
print(tvalues[-1])

gamma = np.array(hydride_interface)
gamma_D1e0 = np.array(hydride_interface_D1e0)
gamma_D1e1 = np.array(hydride_interface_D1e1)
gamma_D1e3 = np.array(hydride_interface_D1e3)
gamma_D1e2 = np.array(hydride_interface_D1e2)


gamma_mu10 = np.array(hydride_interface_mu10)
gamma_mu100 = np.array(hydride_interface_mu100)

c0 = np.array(c0_list)
alpha0 = np.array(alpha0_list)

# Plots for the prediction using constant diffusion coefficient model

# plt.plot(x2_nodes[0:400],alpha_pred[0:400],label='alpha pred', linestyle='dashed',color='red')
# plt.plot(x2_nodes[0:400],alpha[0:400],color='red')

# plt.plot(x2_nodes[0:400],c2_pred[0:400],label='c pred',linestyle='dashed',color='blue')
# plt.plot(x2_nodes[0:400],c2[0:400],color='blue')

# gamma_pred = -1/np.sqrt(3) * np.log(-np.sqrt(3)/(tvalues) * np.log(threshold_alpha))

# for i in range(len(gamma_pred)):
#     if gamma_pred[i] < 0:
#         gamma_pred[i] = 0

# plt.plot(tvalues,hydride_interface,color='black')
# plt.plot(tvalues,gamma_pred,color='black',linestyle='dashed')

print(hydride_interface[-1])
print(tvalues[-1])

# Plots using the new prediction 
 
# plt.plot(tvalues,1/(np.sqrt(3)*(1+np.sqrt(3)*gamma)),linestyle='dashed',color='green',label='c0 pred')
# plt.plot(tvalues,c0,color='green')

# plt.plot(tvalues,alpha0)

tau_c = -np.sqrt(3) * np.log(threshold_alpha)



beta = 1 + tvalues / (np.sqrt(3) * (1-threshold_alpha))

# gamma_pred_leading_order = tvalues / (3*(1-threshold_alpha)) 

# gamma_pred_fake = tvalues / (3*(1-threshold_alpha)) + threshold_alpha * 1/np.sqrt(3) * beta * (np.log(beta) - 1) / (9*np.sqrt(3)*(1-threshold_alpha)**2-np.log(beta))
# gamma_pred = gamma_pred_leading_order  + (threshold_alpha / (9*np.sqrt(3)*(1-threshold_alpha)) * (beta * np.log(beta) - beta) ) / (1-threshold_alpha - threshold_alpha / (9*(1-threshold_alpha)) * np.log(beta))

# gamma_integrate = 1/np.sqrt(3) * np.log(2 * tau_c / (tau_c - np.sqrt(3) * (tvalues - tau_c)))
# x_vals = np.linspace(0,5,100)

# tau_vals_2 = (2*(1-np.exp(-np.sqrt(3)*x_vals))-np.sqrt(3)*x_vals * np.exp(-np.sqrt(3)*x_vals)) * tau_c / np.sqrt(3) + tau_c
# tau_vals = (2-( np.sqrt(3)*x_vals + 2 ) * np.exp(-np.sqrt(3)*x_vals) ) * tau_c / np.sqrt(3)

# print(tau_vals)

# plt.plot(tvalues,gamma)
# plt.plot(tvalues,gamma_pred,linestyle='dashed',label='first correction')


 # GOOD PLOTS

# plt.plot(tvalues,hydride_interface,label='numerics',color='orange')
# plt.plot(tvalues,(-1+np.sqrt(1+2*(tvalues-tau_c)/tau_c))/np.sqrt(3),label='prediction',linestyle='dashed',color='orange')

# plt.xscale('log')
# plt.yscale('log')

# plt.plot(tau_vals_2,x_vals,linestyle='dotted',label='integrated2')

# plt.plot(tvalues,gamma_pred_leading_order,linestyle='dotted',label='leading order')
# plt.plot(tvalues,tvalues / (3*(1-threshold_alpha)))

# plt.xlabel(r'$\tau$')
# plt.ylabel(r'$\Gamma(\tau)$')
# plt.ylabel(r'$c$')

def gam_pred(t):
    if t < tau_c:
        return 0

    return (-1+np.sqrt(1+2*(t-tau_c)/tau_c))/np.sqrt(3)


def alpha_less_than(t,x):
    return np.exp(-tau_c * (1/np.sqrt(3)+(gam_pred(t)-x)))

def alpha_greater_than(t,x):
    return np.exp(-tau_c/np.sqrt(3) * np.exp(np.sqrt(3)*(gam_pred(t)-x)))


def alph_pred(t,x):
    
    if x > gam_pred(t):
        return alpha_greater_than(t,x)
    else:
        return alpha_less_than(t,x)
    
def con_pred(t,x):
    if t < tau_c:
        return 1/np.sqrt(3)
    if x < gam_pred(t):
        return 1/(np.sqrt(3)*(1+np.sqrt(3)*gam_pred(t)))
    else:
        return np.exp(np.sqrt(3)*(gam_pred(t)-x))/(np.sqrt(3)*(1+np.sqrt(3)*gam_pred(t)))

alpha_pred = np.copy(x2_nodes)

for i in range(len(alpha_pred)):
    alpha_pred[i] = alph_pred(1,x2_nodes[i])




# plt.plot(x2_nodes[0:80],alpha[0:80],color='green',label=r'Numerical $\alpha$')
# plt.plot(x2_nodes[0:80],alpha_pred[0:80],linestyle='dashed',color='green',label=r'Predicted $\alpha$')
# plt.plot(hydride_interface[-1] * np.ones(100),np.linspace(0.9,1,100),color='blue',linestyle='solid',label=r'Numerical $\Gamma(\tau=0.875)$')
# plt.plot(gam_pred(t) * np.ones(100),np.linspace(0.9,1,100),color='blue',linestyle='dashed',label=r'Predicted $\Gamma(\tau=0.875)$')

# plt.xlabel(r'z')
# plt.ylabel(r'$\alpha(\tau=0.875,z)$')

conc_pred = np.copy(x2_nodes)

for i in range(len(conc_pred)):
    conc_pred[i] = con_pred(1,x2_nodes[i])

n=350
# plt.plot(x2_nodes[0:n],alpha_pred[0:n],label=r'Predicted concentration',linestyle='dashed',color='red')
# plt.plot(x2_nodes[0:n],interpolated_a_cs_1(x2_nodes)[0:n],label=r'Numerical concentration ($c_s=0.1$)',color='green')
# plt.plot(x2_nodes[0:n],interpolated_a_1(x2_nodes)[0:n],label=r'Numerical concentration ($c_s=0$)',color='blue')

# plt.plot(gam_pred(t) * np.ones(100),np.linspace(0.0,0.1,100),color='blue',linestyle='dashed',label=r'Predicted $\Gamma(\tau)$')

gamma_pred = np.copy(tvalues)
for i in range(len(tvalues)):
    gamma_pred[i] = gam_pred(tvalues[i])

L =  4.381344024813161e-06 * 1e7
T =  56.190976594794634

tvalues = tvalues * T

plt.plot(tvalues,L*gamma,color='magenta',linestyle='dashdot')
plt.plot(tvalues,L*gamma_D1e0,color='red',linestyle='dashed')
plt.plot(tvalues,L*gamma_D1e1,color='red',linestyle='dashed')
plt.plot(tvalues,L*gamma_D1e2,color='red',linestyle='dashed')
plt.plot(tvalues,L*gamma_D1e3,color='red',linestyle='dashed')

plt.plot(tvalues,L*gamma_mu10,color='blue',linestyle='dotted')
plt.plot(tvalues,L*gamma_mu100,color='blue',linestyle='dotted')


plt.plot(tvalues,L*gamma_pred,color='green')




# plt.plot(tvalues,gamma_pred,linestyle='dashed',color='red')

# plt.plot(tvalues,c0)

# plt.xlabel(r'$z$')
# plt.ylabel(r'$c-c_s$')

# plt.plot(tvalues,c0)

# plt.plot(x2_nodes[0:300],alpha[0:300])
# plt.savefig('formal_work/1D_code/regime2a_plots/cs-vary-gamma-plot.eps',format='eps',bbox_inches='tight')
# plt.savefig('formal_work/1D_code/regime2a_plots/pre-fracture.eps',format='eps',bbox_inches='tight')
plt.savefig('formal_work/1D_code/regime2a_plots/figures_latex/D-vary-plot.eps',format='eps',bbox_inches='tight')

# plt.legend()
plt.show()