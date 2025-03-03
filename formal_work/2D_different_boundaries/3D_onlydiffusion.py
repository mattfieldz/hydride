import numpy as np
import scipy as sp
import time

n2 = 401 # metal
m = 401 # width nodes

# number of variables in the bulk
nv = 3  # [H] and [U] in the bulk plus diffusivity of mixture

# dimensional quantities : here approximated for room temp + some ad-hoc choices
L2star = 150 * 1e3 * 1.0e-7  # bulk domain lengthscale 1000um
L3star = 20 * 1.0e-7  # oxide domain lengthscale 20nm
Lmstar = 150 * 1e3 * 1.0e-7 # 500um
D1star = 1.0e-13  # cm^2/s, diffusion of H in UH3  -- *ad hoc* value
D2star = 5.0e-10  # cm^2/s, diffusion of H in U (room temp value)
D3star = 1.0e-13  # cm^2/s, diffusion of H in UO2 (room temp value)
N2star = 8.01e-2  # mol/cm^3, number density of U
Csstar = 1.0e-5  # mol/cm^3, saturation [H] in U -- *ad-hoc* value
Castar = 1.0e-4  # mol/cm^3, surface value for [H] -- *ad-hoc* value

# fixed reference values to keep time/lengthscales consistent in comparisons

Lrefstar = 1e2 * 1.0e-7  # 100nm
# Lrefstar = L2star

Drefstar = D2star  # using the U value

# non-dimensional domains
L3 = L3star / Lrefstar  # oxid
L2 = L2star / Lrefstar  # bulk
Lm = Lmstar / Lrefstar # width

# Choose computation coord Xmax such that x_from_X(Xmax)=L2




# x2_nodes = np.linspace(0,L2,n2)
# x2d_nodes = np.ones(n2)
# Xmax = 1

m2_nodes = np.linspace(0,Lm,m)
m2d_nodes = np.ones(m)

def x_from_X(X):
    BX = 0.01
    return (X + BX) * (X + BX) - BX * BX, 2 * (X + BX)

Xmax = sp.optimize.root(lambda X: x_from_X(X)[0] - L2, 5).x[0]
assert abs(x_from_X(Xmax)[0] - L2) < 1.0e-8
X2_nodes = np.linspace(0, Xmax, n2)  # bulk nodes in comp. coordinate
x2_nodes, x2d_nodes = x_from_X(X2_nodes)





# non-dimensional max = 10^4 sec
tmax = 1e10

np.savetxt(
            f"formal_work/data2D/k0/no_oxide/glascott/x2_nodes_only_7.dat",
            x2_nodes,
        )
# non-dimensional variables
c2 = np.zeros((n2,m), dtype=np.float64)  # metal
alpha = np.ones((n2,m), dtype=np.float64)  # volume fraction of metal
D3 = D3star / Drefstar
D2 = D2star / Drefstar
D1 = D1star / Drefstar
cs = Csstar / Castar  # saturation value \in [0,1)
eps = Castar / N2star  # relative forcing measure
reactK = (
    0  # ad-hoc value, SSI 2024 paper suggests 1.e^4-1.e^5 but based on 1nm scale!
)
D = D2 * np.ones((n2,m), dtype=np.float64)  # mixture diffusion in the bulk, initially =D2

# information output sanity check
print(f"Diffusivity ratio D1*/Dref*={D1star/Drefstar}")
print(f"Diffusivity ratio D2*/Dref*={D2star/Drefstar}")
print(f"Diffusivity ratio D3*/Dref*={D3star/Drefstar}")
print(f"Lengthscale ratio L2*/Lref*={L2star/Lrefstar}")
print(f"Lengthscale ratio L3*/Lref*={L3star/Lrefstar}")
print(f"Diffusion timescale for L3* of oxide T3*={L3star**2/D3star} sec")
print(f"Diffusion timescale for Lref* of metal T2*={Lrefstar**2/D2star} sec")
print(f"Diffusion timescale for Lref* of hydride T1*={Lrefstar**2/D1star} sec")
print(f"Ratio of H to U concentration eps={eps}")
print(f"Relative solubility limit ={cs}")
print(f"Non-dimensional max time ={tmax}")
# metric_file = open("1Dexamples/data/metric.dat", "w")

# step sizes
h2 = L2 / (Xmax*(n2 - 1))
h2s = h2 * h2
hl = Lm / ((m-1))
hls = hl * hl

print(f"grid spacing h2={h2} and hs = {hl}")

dt = 100



# time step the system
t = 0
step = 1  # periodic output counter

computer_t0 = time.time()

def surface_function(m,M):
    radius = 15 * 1e3 * 1e-7
    radius_non_dim = radius / Lrefstar
    index = np.argmin(np.abs(radius_non_dim - m2_nodes))
    return index
radial_value = surface_function(m,Lmstar)
print('radial_value = ',surface_function(m,Lmstar))





while t < tmax:
    # evaluation at the new time step
    t += dt
    # store previous time step solution profiles
    c2_old = np.copy(c2)
    alpha_old = np.copy(alpha)
    D_old = np.copy(D)
    row = []
    col = []
    val = []
    b = []
    count = 0

    for m_i in range(0,m):
        m_s = (n2)*m_i # current node
        m_j = (n2) # spacing
        m_r = m_i * hls * 2 


        for j in range(0, n2):
            k = j + m_s

            if m_i == 0:
                    
                row += 3 * [k]
                col += [(k+2*m_j), (k+m_j), (k)]
                val += [-1, 4, -3]
                b += [0]
                    
                   
            elif m_i == m-1:
                
            
                row += 3 * [k]
                col += [(k-2*m_j), (k-m_j), (k)]
                val += [-1, 4, -3]
                b += [0]
                
            else:    
            
                # at the first and last node we need BCns for the H diffusion
                if (j == 0) or (j == n2 - 1):
                    if j == 0:
                        # at the first bulk node we impose the flux condition
                        if (m_i > 0) and (m_i <= radial_value):
                            row += 3 * [k]
                            col += [
                                k + 2,
                                k + 1,
                                k, 
                            ]
                            val += [-1 * D2 / (2 * h2 * x2d_nodes[0]),   
                                    4 * D2 / (2 * h2 * x2d_nodes[0]),
                                    -3 * D2 / (2 * h2 * x2d_nodes[0])]
                            b += [D3 * (-1) / L3]

                        else:
                            row += 3*[k]
                            col += [k + 2,k + 1,k]
                            val += [-1, 4, -3]
                            b += [0]
                            
                            
                                   
                            
                            
                        
                    if j == n2 - 1:
                        # at the last node we impose no H flux out of the domain
                        row += 3 * [k]
                        col += [(k-2), (k-1), (k)]
                        val += [-1, 4, -3]
                        b += [0]
                else:  # we impose the diffusion equation at internal nodes
                    # these values appear in the diffusion w.r.t. comp. coord.
                    x2dph = 0.5*(x2d_nodes[j+1] + x2d_nodes[j])
                    x2dmh = 0.5*(x2d_nodes[j] + x2d_nodes[j-1])
                    
                    row += 5 * [k]
                    col += [k-1, k, k+1, k - m_j, k + m_j]
                    val += [D2/(h2s*x2d_nodes[j] * x2dmh),-2*D2/hls - 2*D2/(h2s*x2d_nodes[j]) * 0.5*(1/x2dmh + 1/x2dph) - 2/dt, D2/(h2s*x2d_nodes[j]*x2dph), D2/hls - D2/m_r, D2/hls + D2/m_r]
                    b += [- D2/(h2s*x2d_nodes[j]) * (c2_old[j-1,m_i]/x2dmh - 2*c2_old[j,m_i] * 0.5 * (1/x2dph + 1/x2dmh) + c2_old[j+1,m_i]/x2dph)
                          - 2/dt*c2_old[j,m_i] 
                          - D2/hls * (c2_old[j,m_i-1] - 2*c2_old[j,m_i] + c2_old[j,m_i+1])
                          - D2/m_r*(c2_old[j,m_i+1]-c2_old[j,m_i-1])
                    ]
                    
    a = sp.sparse.coo_matrix((val, (row, col)), shape=((n2) * m, (n2) * m))
        # system = petsc.PETScSparseLinearSystem(a, b)
        # x = system.linear_solve()
        # plt.imshow(a.toarray().astype(bool))
        # plt.show()
        # print('det',np.linalg.det(a.toarray()))
    x = sp.sparse.linalg.spsolve(a.tocsr(),b)
    

            
    
    
    
    x_reshaped = x.reshape((n2,m),order='F')
    
    c2 = x_reshaped

    if step % 10 == 0:
        print(
            f"Plot at t={t} tstar={t*(Lrefstar**2)/(Drefstar)} (sec) c_int={c2[0,0]} c_end={c2[-1,0],c2[0,-1]} computer_time ={time.time()-computer_t0}"
        )
        np.savetxt(
            f"formal_work/data2D/k0/no_oxide/glascott/c2_401_only_diff_7{t:.2f}.dat",
            c2,)
    step += 1