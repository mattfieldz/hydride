import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

n_x_nodes = 200



Lstar = 75 * 1e0 * 1e-7 # 1um

D_U_star = 3.0e-19 # guess cm2/s
D_UH3_star = 1e-25 # guess cm2/s

half_life = 12.32 # years
decay_constant = np.log(2) / (half_life*(365 * 24 * 3600)) # decay constant
Nstar = 4.5 * 1e-2 # mol / cm^3 molar density of uranium hydride

Drefstar = D_UH3_star





def compute_released(grams,loading,tritium_percent):
    D_U = D_U_star / Drefstar
    D_UH3 = D_UH3_star / Drefstar

    lambda_ = decay_constant * Lstar**2 / (Drefstar)


    c = np.zeros(n_x_nodes)
    alpha = 0

    D = alpha * D_U + (1-alpha) * D_UH3

    tstar = Lstar ** 2 / Drefstar

    print('lambda : ',lambda_)
    print('Decay constant : ',decay_constant)
    print('D : ', D)
    print('tstar : ',tstar)

    tmax = 16 * 3600 * 24 * 365 / tstar
    t = 0
    step = 0


    x_nodes = np.linspace(0,1,n_x_nodes)
    # h = 1 / n
    h = x_nodes[1]-x_nodes[0]

    h2 = h*h
    dt = 3600 * 8 / tstar

    released = [0]
    generated = [0]
    released_predicted = [0]
    t_days = [0]


    while t < tmax: # loops over time
        row = []
        col = []
        val = []
        b = []

        t += dt
        step += 1
        alpha = (1-np.exp(-lambda_*t))
        G = 3 * lambda_ * (1-alpha)
        

        for j in range(0,n_x_nodes): # loops over space

            if j == 0 or j == n_x_nodes-1: # boundary conditions
                if j==0: # zero flux
                    row += 3 * [j]
                    col += [0,1,2]
                    val += [-3,4,-1]
                    b += [0]
                if j == n_x_nodes-1:
                    row += [j]
                    col += [j]
                    val += [1]
                    b += [0]
                    
            else: # interior nodes
                if x_nodes[j] > 0.0: 
                    row += 3 * [j]
                    col += [j-1,j,j+1]
                    val += [D / h2, -2 * D / h2 - 1/dt, D / h2]
                    b += [-G - c[j]/dt]
                # else: # puts in some region of no source (solid uranium)
                #     row += 3 * [j]
                #     col += [j-1,j,j+1]
                #     val += [D_U / h2, -2 * D_U / h2 - 1/dt, D_U / h2]
                #     b += [- c[j]/dt]
    # Unloaded D
        # D = alpha * D_U + (1-alpha) * D_UH3
    # Loaded D
        D = tritium_percent * (alpha * D_U + (1-alpha) * D_UH3) + (1-tritium_percent) * D_UH3
        
        A = sp.sparse.coo_matrix((val, (row, col)), shape=(n_x_nodes,n_x_nodes))
        
        c_new = sp.sparse.linalg.spsolve(A.tocsr(),b)
        c = c_new
        
        

        if step % 1000 == 0:
            
            interpolated_c = sp.interpolate.CubicSpline(x_nodes, c)
            
            integral_c = interpolated_c.integrate(0.0,1)
            
            c_predicted = np.zeros(len(x_nodes))
            released_predicted_i = 3*(1-np.exp(-lambda_*t))
            for n in range(1000):
                coeff = ((2*n+1)*np.pi/2)
                b_n = 4 * (-1)**(n) / ((2*n+1)*np.pi)
                lambda_n = coeff**2
                c_predicted += 3*b_n * lambda_ * (np.exp(-lambda_n*D_UH3 * t) - np.exp(-lambda_ * t)) / (lambda_ - lambda_n * D_UH3) * np.cos(coeff * x_nodes)
                released_predicted_i -= 24 * lambda_ * (np.exp(-lambda_n*D_UH3 * t) - np.exp(-lambda_ * t)) / ((lambda_ - lambda_n * D_UH3) * (2*n+1)**2*np.pi**2 )
            released_predicted += [released_predicted_i]
            released += [3 * alpha - integral_c]
            generated += [3 * alpha] 

            t_days += [t * tstar/(24*3600)]
            # plt.plot(x_nodes,c)
            # plt.plot(x_nodes,c_predicted,label='predicted')
            # plt.legend()
            # plt.show()
            # np.savetxt(
            #     f"formal_work/data1D/asymptotic/alpha_1001_{t:.5f}.dat",
            #     np.transpose(np.array([t_days, alpha])),
            # )
            
            
            # plt.pause(0.1)
            # plt.plot(c)
            # plt.pause(0.1)

    generated = grams * 2476 * np.array(generated)
    released = grams * 2476 * np.array(released)
    released_predicted = grams * 2476 * np.array(released_predicted)
            
    return t_days, generated, released, released_predicted,

t_days,generated,released,released_predicted = compute_released(0.703,1,0.879)
data_09 = compute_released(0.938,1,0.425)
# data_036 = compute_released()

for i in range(5):
            original_experiments_data = np.loadtxt(f'Helium/diffusion_model/data/original_experiments_{i}.txt')
            plt.scatter(original_experiments_data[0,:],original_experiments_data[1,:],marker='x')

plt.plot(t_days, released,label='0.7 grams / 88 percent',color='blue',linestyle='dashed')
plt.plot(t_days, generated,label='0.7 grams / 88 percent',color='blue',linestyle='dotted')
# plt.plot(t_days,released_predicted,label='predicted release')

plt.plot(data_09[0],data_09[1],label='0.9 grams / 40 percent',linestyle='dotted',color='red')
plt.plot(data_09[0],data_09[2],label='0.9 grams / 40 percent',linestyle='dashed',color='red')


plt.legend()
plt.show()