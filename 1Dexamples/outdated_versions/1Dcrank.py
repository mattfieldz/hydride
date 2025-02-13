
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
# import deps.sp2petsc as petsc

# indexing is 1=Hydride (not present), 2=Metal, 3=Oxide
# discretisation
n3 = 51  # oxide points
n2 = 201 # metal points

# e1 = 4 to jump to concentration nodes
e1 = 2

def x_from_X(X):
    """Defines the (non-uniformly distributed) physical coordinate 'x'
    in terms of the (uniformly distributed) computational coordinate
    'X'. Returns: not only x(X) but also x'(X) and x''(X).

    """
    BX = 0.01  # smaller BX leads to more NON-uniformity
    # +0*X below prevents numpy from returning a scalar when acting
    # on an iterable container for X
    return (X + BX) * (X + BX) - BX * BX, 2 * (X + BX)
    # return X, 1





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
L2 = 1

tmax = 1000

Xmax = sp.optimize.root(lambda X: x_from_X(X)[0] - L2, 5).x[0]
# exit if sanity check fails
assert abs(x_from_X(Xmax)[0] - L2) < 1.0e-8
# nodes in the oxide and bulk
x3_nodes = np.linspace(-L3, 0, n3)  # oxide, the interface is at x=0
X2_nodes = np.linspace(0, Xmax, n2)  # bulk nodes in comp. coordinate

print(Xmax)
# physical node locations
x2_nodes, x2d_nodes = x_from_X(X2_nodes)

# information sanity check
print(f"Diffusivity ratio D3*/D2*={D3star/D2star}")
print(f"Lengthscale ratio L3*/L2*={L3star/L2star}")
print(f"Diffusion time scale for oxide T1*={L3star**2/D3star} sec")
print(f"Diffusion time scale for metal T2*={L2star**2/D2star} sec")

# non-dimensional variables
c3_nodes = np.zeros(n3*e1, dtype=np.float64)
c3_nodes[0] = 1
c2_nodes = np.zeros(n2*e1, dtype=np.float64)
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


Dref = D2star
Lref = L2star

# k_2 = k2star*epsilon**2*N2star**3*Lref**2/Dref
k_2 = 1e6


# analogue to stefan parameter dx/dt = A dc/dx
A = D2star * epsilon / (3)

print('A : ',A)
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

iterations = 0
count_spline = 0


while t < tmax:
    # Tolerance 
    tau = 1
    c3 = np.copy(c3_nodes)
    c2 = np.copy(c2_nodes)

    c3old = np.copy(c3_nodes)
    c2old = np.copy(c2_nodes)

    if iterations > 10:
        dt = 0.005
    iterations = 0

    while tau > 1e-8:  

        


        row = []
        col = []
        val = []
        b = []  # RHS vector for Ax=b

        # surface has C3=1 for j=0
        row += [0,1]
        col += [0,1]
        val += [1,1]
        # b += [1-np.exp(-t)-c3[0],0]
        b += [1 - np.exp(-t) - c3[0],0]
        # surface a1 = 0, a2 = 0, a3 = 1
        





        # internal nodes
        for j in range(1, n3 - 1):
            # concentration nodes
            row += 3 * [e1*j] # this is [j,j,j]
            col += [e1*(j-1), e1*j, e1*(j+1)]
            val += [D3 / h3s, -2 * D3 / h3s - 1 / dt, D3 / h3s]
            b += [-c3_nodes[e1*j] / dt - D3/h3s * c3[e1*(j-1)] + (1/dt+2*D3/h3s) * c3[e1*j]
                  - D3/h3s * c3[e1*(j+1)]]
            # b+=[0]
            # metal nodes
            row += [e1 * j + 1]
            col += [e1 * j + 1]
            val += [1]
            b += [0]

        # interface has equal concentrations: C3 = C2
        j = n3 - 1
        row += 2 * [e1*j]
        col += [e1*j, e1*j + e1]
        val += [1, -1]
        b += [-c3[e1*j] + c2[0]]
    # oxide region
        row += [e1 * j + 1]
        col += [e1 * j + 1]
        val += [1]
        b += [0]

        # interface has matched flux: D3*C3' - D2*C2' = 0 
        j = n3
        D = D2 * (1-c2[1]) + D1 * c2[1]
        
        row += 6 * [e1*j]
        col += [e1*(j - 3),e1*(j - 2), e1*(j - 1), e1*j, e1*(j+1), e1*(j+2)]  # C3_{n3-2},C3_{n3-1},C2_{0},C2_{1}
        # note: only first order accurate, but can use 3 points either side for second order
        val += [D3 / h3 * 0.5,
                 -4 * D3 / h3 * 0.5,
                   3 * D3 / h3 * 0.5,
                 3 * D / (h2*x2d_nodes[0]) * 0.5,
                   -4 * D / (h2*x2d_nodes[0]) * 0.5,
                   D / (h2*x2d_nodes[0]) * 0.5]
        b += [-D3/h3 * 0.5 * (c3[e1*(j-3)] - 4*c3[e1*(j-2)] + 3*c3[e1*(j-1)]) -
               D/(h2*x2d_nodes[0]) * 0.5 * (c2[4] - 4*c2[2] + 3*c2[0])]


        # Alpha values on the interface
        row += 2 * [e1 * j + 1]
        col += [e1 * j,e1 * j + 1]
        val += [3*epsilon * k_2 *(c2[0]-c_sol)**2 * (1-c2[1]),1/dt + epsilon * reac_func(c2[0],k_2,c_sol)]
        b += [(c2_nodes[1]-c2[1])/dt + (1-c2[1]) * epsilon * reac_func(c2[0],k_2,c_sol)]

        # internal nodes in metal region
        for j in range(1, n2 - 1):  # logical node in metal
            k = j + n3  # offset by oxide region eqns

            D_vector = D2 * (1-c2[1::2]) + D1 * c2[1::2]
            D_pos_half = 0.5 * (D_vector[j+1]+D_vector[j])
            D_neg_half = 0.5 * (D_vector[j]+D_vector[j-1])
            x2dph = 0.5 * (x2d_nodes[j + 1] + x2d_nodes[j])
            x2dmh = 0.5 * (x2d_nodes[j] + x2d_nodes[j - 1])

            row += 4 * [e1 * k]
            col += [e1 * (k-1), e1*k,e1*k+1, e1*(k+1)]
            val += [D_neg_half / (h2s * x2d_nodes[j] * x2dmh),
                     -(D_pos_half/x2dph + D_neg_half/x2dmh) / (h2s * x2d_nodes[j]) - 1 / dt - 
                    9 * (c2[e1*j]-c_sol)**2*k_2*(1-c2[e1*j+1]) * np.heaviside(c2[e1*j]-c_sol,0.5),
                    -3*reac_func(c2[j],k_2,c_sol),
                     D_pos_half / (h2s * x2d_nodes[j] * x2dph)]
            b += [-c2_nodes[e1 * j] / dt + 
                  3*reac_func(c2[e1*j],k_2,c_sol)*(1-c2[e1*j+1])
                  -(D_neg_half / (h2s * x2dmh * x2d_nodes[j])) * c2[e1*(j-1)] + 
                  (1/dt + (D_pos_half/x2dph + D_neg_half/x2dmh) / (h2s*x2d_nodes[j]))* c2[e1*j] 
                  - (D_pos_half / (h2s* x2dph * x2d_nodes[j])) * c2[e1*(j+1)]]
            # b += [0]
            # b += [-c2_nodes[e1*j]/dt]

            # Alpha values in internal 
            row += 2*[e1 * k + 1]
            col += [e1 * k ,e1 * k + 1]
            val += [3*epsilon * k_2 *(c2[e1*j]-c_sol)**2 * (1-c2[e1*j+1]),1/dt + epsilon * reac_func(c2[e1*j],k_2,c_sol)]
            b += [(c2_nodes[e1*j+1] - c2[e1*j+1])/dt + (1 - c2[e1*j+1]) * epsilon * reac_func(c2[e1*j],k_2,c_sol)]
                
        # interior boundary has zero first order deriv 
        j = n2 - 1  # logical node in metal
        k = j + n3  # row in mtx
        row += 3 * [e1 * k]
        col += [e1 * (k-2), e1 * (k -1), e1 * (k)]
        val += [-1, 4, -3]
        b += [c2[e1*(j-2)] + 3*c2[e1*j] - 4*c2[e1*(j-1)]]
        


        row += 3*[e1 * k + 1]
        col += [e1 * (k-2) + 1, e1* (k-1) + 1, e1 * (k) + 1]
        val += [-1,4,-3]
        b += [c2[e1*(j-2)+1] + 3*c2[e1*j+1] - 4*c2[e1*(j-1)+1]]

        a = sp.sparse.coo_matrix((val, (row, col)), shape=(e1 * (n3 + n2), e1 * (n3 + n2)))

        # system = petsc.PETScSparseLinearSystem(a, b)

        # x = system.linear_solve()
        
        c_tilde = sp.sparse.linalg.spsolve(a.tocsr(),b)

        
        # Update Conc
        # c2[0::2] += c_tilde[n3*e1:(n2+n3)*e1][::2]
        # c3[0::2] += c_tilde[0:n3*e1][::2]

        # # Update alpha
        # c2[1::2] += c_tilde[n3*e1:(n2+n3)*e1][1::2]
        # c3[1::2] += c_tilde[0:n3*e1][1::2]
        
        c2 += c_tilde[n3*e1:(n2+n3)*e1]
        c3 += c_tilde[0:n3*e1]
        
        tau = np.linalg.norm(c_tilde)
        # print('tau',tau)
        # print('oxide',np.linalg.norm(c_tilde[n3*e1:(n2+n3)*e1][::2]))
        # print('metal',np.linalg.norm(c_tilde[0:n3*e1][::2]))
        iterations += 1
        
        
    # store the solution for this time level
    
    c3_nodes[0:n3*e1] = c3
    c2_nodes[0:n2*e1] = c2
    x = np.zeros((n2+n3)*e1)
    x[0:n3*e1] = c3
    x[n3*e1:(n2+n3)*e1] = c2
    coords = np.zeros((n2+n3))
    coords[0:n3] = x3_nodes
    coords[n3:(n2+n3)] = x2_nodes

    t += dt
    tspace.append(t*Lref**2/(Dref))
    interfaceconc.append(c2_nodes[0])
    interfacealpha.append(1-c2_nodes[1])
    # save every 10 steps

    if step % 100 == 0:
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
        # while alph > 0.02:
        #     alph = c2[e1*index+1]
        #     index += 1
        #     if e1*index == len(c2):
        #         index = 0
        #         break 
        
        # x1.append(coords[n3+index-1])
        
        # computing the maximum reaction term position
        spline = sp.interpolate.InterpolatedUnivariateSpline(x2_nodes,c2[0::2]**3 * (1-c2[1::2]),k=4)
        
        

        cr_pointsh = spline.derivative().roots()
    
        cr_valuesh = spline(np.append(cr_pointsh,(0,L2)))
        cr_maximal_indexh = np.argmax(cr_valuesh)
        x1_value = np.append(cr_pointsh,(0,L2))[cr_maximal_indexh]

        if len(spline.derivative().roots()) == 1:
            count_spline += 1

        if x1_value > 0.99:
            x1.append(0)
        else:
            x1.append(x1_value)
        print(x1_value)

        axis[0,0].cla()
        axis[0,1].cla()
        axis[1,0].cla()
        axis[1,1].cla()

        axis[0,0].set_title('Hydride thickness')
        axis[0,1].set_title('UH3 concentration')
        axis[1,0].set_title('Reaction term')
        axis[1,1].set_title('UH3 conc @ interface')
        
        x1_100 = (np.array(x1)*Lref)
        x1_array_log = np.log10(np.array(x1)*Lref)
        log_tspace_100 = np.log10(np.array(tspace)[::100])
        tspace_100 = np.array(tspace)[::100]
        pcm_conc = axis[0,0].plot(log_tspace_100,x1_array_log)
        plot_A = axis[0,0].plot(log_tspace_100,np.log10(np.sqrt(2*A*(1-c_sol)))+0.5*log_tspace_100)

        # pcm_conc = axis[0,0].plot(tspace_100,x1_100)
        # plot_A = axis[0,0].plot(tspace_100,np.sqrt(2*A*(1-c_sol)*tspace_100))


        pcm_alpha = axis[0,1].plot(coords,1-x[1::2])
        pcm_beta = axis[1,0].plot(tspace,interfaceconc)
        # pcm_beta = axis[1,0].plot(c2[0::2]**3 * (1-c2[1::2]))
        # pcm_gamma = axis[1,1].plot(tspace,interfacealpha)
        
        plot_reac = axis[1,1].plot(x2_nodes,c2[0::2]**3 * (1-c2[1::2]))


   
        fig.suptitle(t*Lref**2/(Dref))
        plt.pause(0.01)

    step += 1

# np.savetxt(
#             f"1Dexamples/data/fullnumuh3.dat",
#             interfacealpha)
# np.savetxt(
#             f"1Dexamples/data/fullnumtimes.dat",
#             tspace)
# np.savetxt(
#             f"1Dexamples/data/fullnumH.dat",
#             interfaceconc)