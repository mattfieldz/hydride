
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import time
# import deps.sp2petsc as petsc
computer_time0 = time.time()
# indexing is 1=Hydride (not present), 2=Metal, 3=Oxide
# discretisation
n3 = 51  # oxide points
n2 = 101 # metal points

m = 20 # width nodes

# e1 = 4 to jump to concentration nodes
e1 = 2

# dimensional quantities, natchiar 2024 for diffusivity of hydrogen in U, UH3, UO2
L3star = 10*1.0e-7  # oxide lengthscale 10 nm = 1.e-8 m = 1.e-6 cm
L2star = 1000*1.0e-7  # bulk lengthscale  50 um = 5.e-5 m = 5.e-3 cm
Lmstar = 1e6 * 1.0e-7 # millimetre width


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
Lm = Lmstar / L2star
L2 = 1
x3_nodes = np.linspace(0, L3, n3)
x2_nodes = np.linspace(L3, L3 + L2, n2)
tmax = 200

# information sanity check
print(f"Diffusivity ratio D3*/D2*={D3star/D2star}")
print(f"Lengthscale ratio L3*/L2*={L3star/L2star}")
print(f"Diffusion time scale for oxide T1*={L3star**2/D3star} sec")
print(f"Diffusion time scale for metal T2*={L2star**2/D2star} sec")

# non-dimensional variables
c3_nodes = np.zeros((n3*e1,m), dtype=np.float64)
c3_nodes[0] = 1
c2_nodes = np.zeros((n2*e1,m), dtype=np.float64)
c = np.zeros(((n2+n3)*e1,m))
# c2_nodes[1::2] = 1s

D3 = D3star / D2star 
D1 = D1star / D2star 
D2 = 1

# reaction rate

# step sizes
h3 = L3 / (n3 - 1)
h3s = h3 * h3
h2 = L2 / (n2 - 1)
h2s = h2 * h2
hl = Lm / (m-1)
hls = hl * hl

Dref = D2star
Lref = L2star

# k_2 = k2star*epsilon**2*N2star**3*Lref**2/Dref
k_2 = 1e5


print('k_2',k_2)
# dt = Dref/(L2star**2) * epsilon * 10000
dt = 0.02




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
    return k*(conc-c_crit)**3*np.heaviside(conc-c_crit,0.5)
    
tspace = []
alpha_at_boundary = []
x1 = []

t_values_for_plotting = []

while t < tmax:
    # Tolerance 
    tau = 1
    c3 = np.copy(c3_nodes)
    c2 = np.copy(c2_nodes)
    iterations = 0
    t += dt
    
    while tau > 1e-8:  

        


        row = []
        col = []
        val = []
        b = []  # RHS vector for Ax=b

       
        


        for m_i in range(0,m):
            m_s = e1*(n2+n3)*m_i
            m_j = (n2+n3)

            
          
            if (m_i == 0) or (m_i == m-1):
                if m_i == m-1:
                    for j in range(0,n2+n3):
                        row += 3 * [e1 * j + m_s]
                        col += [e1 * (j-2*m_j) + m_s, e1 * (j-m_j) + m_s, e1 * (j)+m_s]
                        val += [-1, 4, -3]
                        b += [c[e1*(j),m_i-2] + 3*c[e1*j,m_i] - 4*c[e1*(j),m_i-1]]
                    
                        row += 3*[e1 * j + 1+m_s]
                        col += [e1 * (j-2*m_j) + 1 + m_s, e1* (j-m_j) + 1 + m_s, e1 * (j) + 1 + m_s]
                        val += [-1,4,-3]
                        b += [c[e1*(j)+1,m_i-2] + 3*c[e1*j+1,m_i] - 4*c[e1*(j)+1,m_i-1]]
                
                if m_i == 0:
                    for j in range(0,n2+n3):
                        row += 3 * [e1 * j + m_s]
                        col += [e1 * (j+2*m_j) + m_s, e1 * (j+m_j) + m_s, e1 * (j)+m_s]
                        val += [-1, 4, -3]
                        b += [c[e1*(j),m_i+2] + 3*c[e1*j,m_i] - 4*c[e1*(j),m_i+1]]
                    
                        row += 3*[e1 * j + 1+m_s]
                        col += [e1 * (j+2*m_j) + 1 + m_s, e1* (j+m_j) + 1 + m_s, e1 * (j) + 1 + m_s]
                        val += [-1,4,-3]
                        b += [c[e1*(j)+1,m_i+2] + 3*c[e1*j+1,m_i] - 4*c[e1*(j)+1,m_i+1]]
            # internal nodes
            else:
                # print(m_s)
                # surface has C3=1 for j=0
                
                row += [0+m_s,1+m_s]
                col += [0+m_s,1+m_s]
                val += [1,1]
                # b += [1-np.exp(-t)-c3[0],0]
                b += [np.exp(-(100*(m_i/m))**2)-c3[0,m_i],0]
                # b += [1-c3[0,m_i],0]
                
                # surface a1 = 0, a2 = 0, a3 = 1

                for j in range(1, n3 - 1):

                

                   
            
                    row += 5 * [e1*j+m_s] # this is [j,j,j]
                    col += [e1*(j-1)+m_s, e1*j+m_s, e1*(j+1)+m_s, e1*(j-m_j) + m_s, e1*(j+m_j) + m_s]
                    val += [D3 / h3s, -2 * D3 / h3s - 2 * D3 / hls - 1 / dt, D3 / h3s, D3 / hls, D3 / hls]
                    b += [-c3_nodes[e1*j,m_i] / dt - D3/h3s * c3[e1*(j-1),m_i] + (1/dt+2*D3/h3s + 2*D3/hls) * c3[e1*j,m_i]
                        - D3/h3s * c3[e1*(j+1),m_i] - D3/hls * c3[e1*(j),m_i-1] - D3/hls * c3[e1*(j),m_i+1]]
                    # b+=[0]
                    # metal nodes
                    row += [e1 * j + 1+m_s]
                    col += [e1 * j + 1+m_s]
                    val += [1]
                    b += [0]



                # interface has equal concentrations: C3 = C2
                j = n3 - 1
                row += 2 * [e1*j+m_s]
                col += [e1*j+m_s, e1*(j+1)+m_s]
                val += [1, -1]
                b += [-c3[e1*j,m_i] + c2[0,m_i]]
            # oxide region
                row += [e1 * j + 1 + m_s]
                col += [e1 * j + 1 + m_s]
                val += [1]
                b += [0]

                # interface has matched flux: D3*C3' - D2*C2' = 0 
                j = n3
                D = D2 * (1-c2[1,m_i]) + D1 * c2[1,m_i]

                row += 6 * [e1*j+m_s]
                col += [e1*(j - 3)+m_s,e1*(j - 2)+m_s, e1*(j - 1)+m_s, e1*j+m_s, e1*(j+1)+m_s, e1*(j+2)+m_s]  # C3_{n3-2},C3_{n3-1},C2_{0},C2_{1}
                # note: only first order accurate, but can use 3 points either side for second order
                val += [D3 / h3 * 0.5, -4 * D3 / h3 * 0.5, 3 * D3 / h3 * 0.5, 3 * D / h2 * 0.5, -4 * D / h2 * 0.5, D / h2 * 0.5]
                b += [-D3/h3 * 0.5 * (c3[e1*(j-3),m_i] - 4*c3[e1*(j-2),m_i] + 3*c3[e1*(j-1),m_i]) -
                    D/h2 * 0.5 * (c2[4,m_i] - 4*c2[2,m_i] + 3*c2[0,m_i])]


                # Alpha values on the interface
                row += 2 * [e1 * j + 1+m_s]
                col += [e1 * j + m_s,e1 * j + 1+m_s]
                val += [3*epsilon * k_2 *(c2[0,m_i]-c_sol)**2 * (1-c2[1,m_i]),1/dt + epsilon * reac_func(c2[0,m_i],k_2,c_sol)]
                b += [(c2_nodes[1,m_i]-c2[1,m_i])/dt + (1-c2[1,m_i]) * epsilon * reac_func(c2[0,m_i],k_2,c_sol)]

                # internal nodes in metal region
                for j in range(1, n2 - 1):  # logical node in metal
                    k = j + n3  # offset by oxide region eqns

                    D_vector = D2 * (1-c2[1::2,m_i]) + D1 * c2[1::2,m_i]
                    D_pos_half = 0.5 * (D_vector[j+1]+D_vector[j])
                    D_neg_half = 0.5 * (D_vector[j]+D_vector[j-1])
                    
                    D_left_right_vector = D2 * (1-c2[j,:]) + D1 * c2[j,:]
                    D_right_half = 0.5 * (D_left_right_vector[m_i+1]+D_left_right_vector[m_i])
                    D_left_half = 0.5 * (D_left_right_vector[m_i]+D_left_right_vector[m_i-1])


                    row += 6 * [e1 * k + m_s]
                    col += [e1 * (k-1) + m_s, e1 * k + m_s, e1* k +1 + m_s, e1* (k+1) + m_s, e1 * (k-m_j) + m_s, e1 * (k+m_j) + m_s ]
                    val += [D_neg_half / h2s,
                            -(D_pos_half + D_neg_half)/h2s - (D_left_half + D_right_half) / hls - 1 / dt - 
                            9 * (c2[e1*j,m_i]-c_sol)**2*k_2*(1-c2[e1*j+1,m_i]) * np.heaviside(c2[e1*j,m_i]-c_sol,0.5),
                            -3*reac_func(c2[j,m_i],k_2,c_sol),
                            D_pos_half / h2s, D_left_half / hls , D_right_half / hls]
                    b += [-c2_nodes[e1 * j,m_i] / dt + 
                        3*reac_func(c2[e1*j,m_i],k_2,c_sol)*(1-c2[e1*j+1,m_i])
                        -(D_neg_half / h2s) * c2[e1*(j-1),m_i] + (1/dt + (D_pos_half + D_neg_half)/h2s + (D_left_half + D_right_half) / hls)* c2[e1*j,m_i] 
                        - (D_pos_half / h2s) * c2[e1*(j+1),m_i]
                        - (D_left_half / hls) * c2[e1*(j),m_i-1] - (D_right_half / hls) * c2[e1*(j),m_i+1] ]
                    # b += [0]
                    # b += [-c2_nodes[e1*j]/dt]

                    # Alpha values in internal 
                    row += 2*[e1 * k + 1+m_s]
                    col += [e1 * k + m_s , e1 * k + 1+m_s]
                    val += [3*epsilon * k_2 *(c2[e1*j,m_i]-c_sol)**2 * (1-c2[e1*j+1,m_i]),1/dt + epsilon * reac_func(c2[e1*j,m_i],k_2,c_sol)]
                    b += [(c2_nodes[e1*j+1,m_i] - c2[e1*j+1,m_i])/dt + (1 - c2[e1*j+1,m_i]) * epsilon * reac_func(c2[e1*j,m_i],k_2,c_sol)]
                        
                # interior boundary has zero first order deriv 
                j = n2 - 1  # logical node in metal
                k = j + n3  # row in mtx
                row += 3 * [e1 * k + m_s]
                col += [e1 * (k-2) + m_s, e1 * (k -1) + m_s, e1 * (k)+m_s]
                val += [-1, 4, -3]
                b += [c2[e1*(j-2),m_i] + 3*c2[e1*j,m_i] - 4*c2[e1*(j-1),m_i]]
                


                row += 3*[e1 * k + 1+m_s]
                col += [e1 * (k-2) + 1 + m_s, e1* (k-1) + 1 + m_s, e1 * (k) + 1 + m_s]
                val += [-1,4,-3]
                b += [c2[e1*(j-2)+1,m_i] + 3*c2[e1*j+1,m_i] - 4*c2[e1*(j-1)+1,m_i]]

        
        
        # print(len(row),len(col),len(val),e1*(n3*n2)*m)
        # print(max(row))



        a = sp.sparse.coo_matrix((val, (row, col)), shape=(e1 * (n3 + n2) * m, e1 * (n3 + n2)*m))
        
        c_solve = sp.sparse.linalg.spsolve(a.tocsr(),b)
        
        c_tilde = c_solve.reshape((e1*(n2+n3),m),order='F')

        c2 += c_tilde[n3*e1:(n2+n3)*e1,:]
        c3 += c_tilde[0:n3*e1,:]
        
        c[0:n3*e1,:] = c3
        c[n3*e1:(n2+n3)*e1,:] = c2

        tau = np.linalg.norm(c_tilde)
  
        iterations += 1 
        
    # store the solution for this time level
    
    c3_nodes[0:n3*e1,:] = c3
    c2_nodes[0:n2*e1,:] = c2
    x = np.zeros(((n2+n3)*e1,m))
    x[0:n3*e1] = c3
    x[n3*e1:(n2+n3)*e1] = c2
    
    
    tspace.append(t*Lref**2/(Dref))
    interfaceconc.append(c2_nodes[0,0])
    interfacealpha.append(1-c2_nodes[1,0])
    
    # save every 10 steps

    if step % 100 == 0:
        print('time to compute',time.time()-computer_time0)
        print('iterations',iterations)
        print(
            f"Save at t={t} c_int={c3_nodes[e1*(n3-1)]}={c2_nodes[0]} c_end={c2_nodes[n2-1]}"
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

        # alph = 1
        # index=0
        # while alph > 0.01:
        #     alph = c2[e1*index+1]
        #     index += 1
        #     if e1*index == len(c2):
        #         index = 0
        #         break 
        
        
            
        
        t_values_for_plotting.append(t)

        length_vector = np.zeros(n2+n3)
        length_vector[0:n3] = np.linspace(-L3,0,n3)
        length_vector[n3:n3+n2] = np.linspace(0,L2,n2)
        
        # x1.append(length_vector[n3+index-1])
        # plt.cla()
        # axis[0,0].cla()
        # axis[0,1].cla()
        # axis[1,0].cla()
        # axis[1,1].cla()

        # axis[0,0].set_title('H concentration cross section')
        # axis[0,1].set_title('UH3 concentration cross section')
        # axis[1,0].set_title('H conc ')
        # axis[1,1].set_title('UH3 conc')

        
        # pcm_conc = axis[0,0].plot(length_vector,x[0::2,int(m/2)])
        # pcm_conc = axis[0,0].plot(np.sqrt(np.array(tspace)[::100]),x1)

        # pcm_alpha = axis[0,1].plot(length_vector,1-x[1::2,int(m/2)])
        
        # pcm_beta = axis[1,0].imshow(c2[0::2],aspect='auto',extent=[-Lm/2,Lm/2,-L2,0])
        # pcm_gamma = axis[1,1].imshow(c2[1::2],aspect='auto',extent=[-Lm/2,Lm/2,-L2,0])
        # fig.suptitle(t*Lref**2/(Dref))
        # cb1 = plt.colorbar(pcm_beta,ax=axis[1,0])
        # cb2 = plt.colorbar(pcm_gamma,ax=axis[1,1])
        # plt.pause(0.01)
        # cb1.remove()
        # cb2.remove()
        
        np.savetxt(f"formal_work/data2D/k1e5/c3_t{t:.2f}.dat",c3
        )
        np.savetxt(
            f"formal_work/data2D/k1e5/c2_t{t:.2f}.dat",c2
        )


    step += 1

np.savetxt(
    f"formal_work/data2D/k1e5/tspace.dat",
    np.array(t_values_for_plotting)
)

# np.savetxt(
#             f"1Dexamples/data/fullnumuh3.dat",
#             interfacealpha)
# np.savetxt(
#             f"1Dexamples/data/fullnumtimes.dat",
#             tspace)
# np.savetxt(
#             f"1Dexamples/data/fullnumH.dat",
#             interfaceconc)