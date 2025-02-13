
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
# import deps.sp2petsc as petsc

# indexing is 1=Hydride (not present), 2=Metal, 3=Oxide
# discretisation
  # oxide points
n2 = 1001 # metal points

# e1 = 4 to jump to concentration nodes
e1 = 2

# dimensional quantities, natchiar 2024 for diffusivity of hydrogen in U, UH3, UO2
L3star = 10*1.0e-7  # oxide lengthscale 10 nm = 1.e-8 m = 1.e-6 cm
L2star = 1000*1.0e-7  # bulk lengthscale  50 um = 5.e-5 m = 5.e-3 cm

D3star = 1.18e-12  # cm^2/s, diffusion of H in UO2 @ room temp
D2star = 1.49e-10  # cm^2/s, diffusion of H in U @ room temp
D1star = 1.0e-13


N1star = 4.54e-2 
N2star = 8.01e-2
N3star = 4.06e-2

k2star = 2.0e13

Castar = 1.0e-4
Csstar = 1.0e-5

epsilon = Castar / N2star


c_sol = 0.1

# non-dimensional domains
L3 = L3star / L2star
L2 = 100

x2_nodes = np.linspace(0,L2, n2)
tmax = 2.0

# information sanity check
print(f"Diffusivity ratio D3*/D2*={D3star/D2star}")
print(f"Lengthscale ratio L3*/L2*={L3star/L2star}")
print(f"Diffusion time scale for oxide T1*={L3star**2/D3star} sec")
print(f"Diffusion time scale for metal T2*={L2star**2/D2star} sec")

# non-dimensional variables
c2_nodes = 0*np.ones(n2*e1, dtype=np.float64)
c2_nodes[1::2] = 1

D3 = D3star / D2star
D1 = D1star / D2star
D2 = 1

# reaction rate

# step sizes


h2 = L2 / (n2 - 1)
h2s = h2 * h2

delta_oxide = L3

Dref = D2star
Lref = L2star

# k_2 = k2star*epsilon**2*N2star**3*Lref**2/Dref
k_2 = 1e5


print('k_2',k_2)
print('delta_oxide', delta_oxide)
# dt = Dref/(L2star**2) * epsilon * 10000
dt = 0.001




# time step
t = 0
step = 1

# Plotting
fig,axis = plt.subplots(2,2)
axis[0,0].set_title('H concentration')
axis[0,1].set_title('UH3 concentration')

interfaceconc = []
interfacealpha = []


def reac_func(conc,k,c_crit):
    return (conc)**3*np.heaviside(conc-c_crit,0.5)
    
tspace = []
alpha_at_boundary = []



iterates = 0
while t < tmax:
    iterates = 0
    # Tolerance 
    tau = 1
    c2_update_nodes = np.copy(c2_nodes)
    iterations = 0

    while tau > 1e-10:  

        c2 = c2_update_nodes


        row = []
        col = []
        val = []
        b = []  # RHS vector for Ax=b

        # interface has matched flux: D3*C3' - D2*C2' = 0 
        j = 0
        k=0
        D = D2 * (1-c2[1]) + D1 * c2[1]
        
        row += 4 * [e1*j]
        col += [e1*(j), e1 * (j) + 1,e1*(j + 1), e1*(j + 2)]  # C3_{n3-2},C3_{n3-1},C2_{0},C2_{1}
        # note: only first order accurate, but can use 3 points either side for second order

        deriv_cg = 0.5 / h2 * (-3 * c2[0] + 4 * c2[2] - c2[4])

        val += [-3 * D2 * c2[1] / h2 * 0.5, 
                 D2 * deriv_cg,
                 4 * D2 * c2[1] / h2 * 0.5,
                 -D2 * c2[1] / h2 * 0.5]
        b += [D3 * (c_sol - 1) / delta_oxide - D2 * c2[1] * deriv_cg]

        row += 2 * [e1 * k + 1]
        col += [e1 * k,e1 * k + 1]
        val += [3 * c2[e1*j]**2 * c2[e1*j+1],1/dt + c2[e1*j]**3 ]
        b += [(c2_nodes[e1*j+1] - c2[e1*j+1])/dt - (c2[e1*j+1]) * c2[e1*j]**3]
                


        # Alpha values on the interface
        
        # internal nodes in metal region
        for j in range(1, n2 - 1):  # logical node in metal
            k = j  # offset by oxide region eqns

            D_vector = D2 * (c2[1::2]) + D1 * (1-c2[1::2])
            D_pos_half = 0.5 * (D_vector[j+1] + D_vector[j])
            D_neg_half = 0.5 * (D_vector[j] + D_vector[j-1])
            

            row += 4 * [e1 * k]
            col += [e1 * (k-1), e1*k,e1*k+1, e1*(k+1)]
            val += [D_neg_half / h2s, -(D_pos_half + D_neg_half) / h2s  - 
                    9 * (c2[e1*j])**2*(c2[e1*j+1]), 3*c2[e1*j]**3,
                     D_pos_half / h2s]
            b += [
                  3*c2[e1*j]**3*(c2[e1*j+1])
                  -(D_neg_half / h2s) * c2[e1*(j-1)] + ((D_pos_half + D_neg_half) / h2s)* c2[e1*j] 
                  - (D_pos_half / h2s) * c2[e1*(j+1)]]
            # b += [0]
            # b += [-c2_nodes[e1*j]/dt]

            # Alpha values in internal 
            row += 2 * [e1 * k + 1]
            col += [e1 * k,e1 * k + 1]
            val += [3 * c2_update_nodes[e1*k]**2 * c2_update_nodes[e1*k+1],1/dt + c2_update_nodes[e1*j]**3 ]
            b += [(c2_nodes[e1*j+1] - c2_update_nodes[e1*j+1])/dt - (c2_update_nodes[e1*j+1]) * c2_update_nodes[e1*j]**3]
                
        # interior boundary has zero first order deriv 
        j = n2 - 1  # logical node in metal
        k = j  # row in mtx
        row += 3 * [e1 * k]
        col += [e1 * (k-2), e1 * (k -1), e1 * (k)]
        val += [-1, 4, -3]
        b += [c2_update_nodes[e1*(j-2)] + 3*c2_update_nodes[e1*j] - 4*c2_update_nodes[e1*(j-1)]]
        


        row += 3*[e1 * k + 1]
        col += [e1 * (k-2) + 1, e1 * (k-1)+1, e1 * (k) + 1]
        val += [-1,4,-3]
        b += [c2_update_nodes[e1*(j-2)+1] + 3*c2_update_nodes[e1*j+1] - 4*c2_update_nodes[e1*(j-1)+1]]

        a = sp.sparse.coo_matrix((val, (row, col)), shape=(e1 * (n2), e1 * (n2)))

        # system = petsc.PETScSparseLinearSystem(a, b)

        # x = system.linear_solve()
        # if t == 0:
        #     array = a.toarray()
        #     print(array.astype('i'))
        #     print(np.linalg.det(array))
        #     print(iterates)


        c_tilde = sp.sparse.linalg.spsolve(a.tocsr(),b)
        
        
        # Update Conc
        # c2_update_nodes[0::2] += c_tilde[n3*e1:(n2+n3)*e1][::2]
        # c3_update_nodes[0::2] += c_tilde[0:n3*e1][::2]

        # # Update alpha
        # c2_update_nodes[1::2] += c_tilde[n3*e1:(n2+n3)*e1][1::2]
        # c3_update_nodes[1::2] += c_tilde[0:n3*e1][1::2]
        
        c2_update_nodes += c_tilde
        
        tau = np.linalg.norm(c_tilde)
        # print('tau',tau)
        # print('oxide',np.linalg.norm(c_tilde[n3*e1:(n2+n3)*e1][::2]))
        # print('metal',np.linalg.norm(c_tilde[0:n3*e1][::2]))
        iterates += 1
        
        
    # store the solution for this time level
    # print(iterates)
    
    c2_nodes = np.copy(c2_update_nodes)
    x = np.zeros((n2)*e1)
    
    x = c2_update_nodes
    
    tspace.append(t*Lref**2/(Dref) / (epsilon * k_2**(0.25)))
    interfaceconc.append(c2_nodes[0]*k_2**(-0.25)+c_sol)
    interfacealpha.append(c2_nodes[1])
    # save every 10 steps
    t += dt
    if step % 10 == 0:
        print('iterations',iterations)
        print(
            f"Save at t={t} c_int={c2_nodes[0]} c_end={c2_nodes[n2-1]}"
        )
        print(
            f"percent UH3={np.sum(x[1::2])/n2}"
        )
        # np.savetxt(
        #     f"./c3_{t:.2f}.dat",
        #     np.transpose(np.array([x3_nodes, c3_nodes])),
        # )
        # np.savetxt(
        #     f"./c2_{t:.2f}.dat",
        #     np.transpose(np.array([x2_nodes, c2_nodes])),
        # )
        # axis[0,0].cla()
        # axis[0,1].cla()
        
        
        length_vector= np.linspace(0,L2,n2)
        
        axis[0,0].cla()
        axis[0,1].cla()
        axis[1,0].cla()
        axis[1,1].cla()

        axis[0,0].set_title('H concentration')
        axis[0,1].set_title('UH3 concentration')
        axis[1,0].set_title('H conc @ interface')
        axis[1,1].set_title('UH3 conc @ interface')


        pcm_conc = axis[0,0].plot(length_vector,x[0::2]*k_2**(-0.25)+c_sol)
        pcm_alpha = axis[0,1].plot(length_vector,1-x[1::2])
        pcm_beta = axis[1,0].plot(tspace,interfaceconc)
        pcm_gamma = axis[1,1].plot(tspace,interfacealpha)
        # print(x[0::2])
    # axis.set_title(counter)
        fig.suptitle(t*Lref**2/(Dref) / (epsilon * k_2**(0.25)))
        plt.pause(0.01) 

    step += 1

np.savetxt(
            f"1Dexamples/data/asymptuh3.dat",
            interfacealpha)
np.savetxt(
            f"1Dexamples/data/asympttime.dat",
            tspace)
np.savetxt(
            f"1Dexamples/data/asymptH.dat",
            interfaceconc)