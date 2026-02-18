import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
plt.style.use('formal_work/minorticks.mplstyle')
plt.margins(0)

mu = 33.6
n_mh = 0.06098
b = 2.328
sigma_f = 2.4

gamma = 14.7

f_p = 0.64
R0 = 19.5
V_c = 0.15





def func_v_a(p):
    A = 18.573
    B = -7.100
    C = 5.375
    return A * p**(-1/3) + B * p**(-2/3) + C * p**(-1)



def func_p(r,R):
    alpha = r/(R)
    
    eps = (alpha/(2-alpha))**3
    beta = 3/4 * alpha**2
    beta = 0 * alpha


# Regular pressure
    pressure = 2 * gamma / r + mu * b / (r * (1+eps))


# Modified pressure to include extra effects
    # pressure = np.copy(r)
    # for i in range(len(r)):
    #     if r[i] < 4 * np.pi*b:
    #         pressure[i] = 2 * gamma / r[i] + mu * b / (r[i] * (1+eps[i]) * (1-beta[i]))
    #     else:
    #         pressure[i] = 2 * gamma / r[i] + mu / (4*np.pi)



    return pressure

def func_p_crit(r,R):
    n_b = f_p / (4 * R**3 * np.pi / 3 )
    return 2 * gamma / r + sigma_f * (1/(np.pi * r**2 * n_b**(2/3)) - 1)


# rr = np.linspace(1,100,1000)
# RR = np.linspace(1,500,1000)

# rrr, RRR = np.meshgrid(rr,RR)

# ppp = func_p(rrr,RRR)
# pppccc = func_p_crit(rrr,RRR)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# # ax.plot_wireframe(rrr, RRR, ppp, color='orange')
# ax.plot_wireframe(rrr,RRR,pppccc,color='blue')
# ax.set_xlabel('r')
# ax.set_ylabel('R')

# plt.show()

def direct_comp(r,R):

    # n_b_no_r = f_p / (4/3*np.pi)
    # beta = 1 / (n_b_no_r**(2/3) * np.pi)
    # beta = (3*f_p / 4)**(-2/3) * np.pi**(-1/3)
    # return alpha**5 - 6 * alpha**4 + (12-beta) * alpha**3 + (6*beta - 8) * alpha**2 - (8*mu*b/(R*sigma_f) + 12*beta) * alpha + 8*beta 
    # return mu * b * r**4 - 6 * R * mu * b * r**3 + (12*R**2 * mu * b - 8 * sigma_f * R**3) * r**2 - (8 * R**3 * mu * b) * r + beta * sigma_f * R**5

    return func_p(r,R) - func_p_crit(r,R)

def poly_for_alpha(r,R):
    alpha = r/R
    A = (3 * f_p / (4*np.pi))**(2/3) * np.pi
    B = mu * b /(sigma_f*R)

    # return alpha**4 * (6*sigma_f - B) + alpha**3 * (-12*sigma_f + 6*B) + alpha**2 * (8*sigma_f - 6*sigma_f/A - 12*B) + alpha * (8*B+12*sigma_f/A) + (-8*sigma_f/A)
    return alpha**4 * (6 - B) + alpha**3 * (-12 + 6*B) + alpha**2 * (8 - 6/A - 12*B) + alpha * (8*B+12/A) + (-8/A)

R_test = 40
rr = np.linspace(20,300,100)


# plt.plot(rr,func_p(rr,R_test),label='p')
# plt.plot(rr,func_p_crit(rr,R_test),label='p_c')
# plt.plot(R_test*np.ones(100),np.linspace(-3,3,100),label='R')

# plt.plot(rr,func_p(rr,R_test)-func_p_crit(rr,R_test),label='crit')
# plt.plot(rr,poly_for_r(rr,R_test))


# plt.plot(rr,polynomial_for_r(rr,R_test))
# plt.legend()
# plt.show()


def root_finder():

    n=1000

    cum_diff = 0


    R_array = np.linspace(0.1,10,n)
    r_solution_list = []
    r_other_branch_list = []
    for R_index in range(n):
        R = R_array[R_index]

        r_array = np.linspace(0.001,R+100,n)

        poly_array = R * poly_for_alpha(r_array,R)
        poly_array_test = direct_comp(r_array,R)
        interpolated_r = sp.interpolate.CubicSpline(r_array, poly_array)
        roots = interpolated_r.roots(extrapolate=False)

        interpolated_r_test = sp.interpolate.CubicSpline(r_array, poly_array_test)
        roots_test = interpolated_r_test.roots(extrapolate=False)
        cum_diff += abs(roots_test[0] - roots[0]) / roots_test[0]
        # print(roots)
        # print(R_index)
        print(R,roots/R)
        
        # print(R)
        r_solution_list.append(roots[0])
    # alpha_solution = roots[0]
    print('cum_diff',cum_diff)
    # print(roots)

    # alpha_solution = roots[0]


    r_solution = np.array(r_solution_list)
    # r_other_branch = np.array(r_other_branch_list)
    return R_array, r_solution

def He_M_func(r,R):
    pressure = func_p(r,R)
    pressure_critical = func_p_crit(r,R)
    v_a = func_v_a(pressure)
    He_M = f_p / (n_mh * v_a) * (r/R)**3
    return He_M
def R_dist_func(R):
    
    m = 11 
    sigma = 1.14 * 1
    Rd = R - R0
    return ((2*np.pi)**0.5 * Rd * sigma) ** (-1) * np.exp( -(np.log(Rd/m))**2 / (2*sigma**2))


# R_array, r_solution = root_finder()


# He_M = He_M_func(r_solution,R_array)
# pressure = func_p(r_solution,R_array)
# v_a = func_v_a(pressure)

def plot_polynomial():
    sing = mu * b /(sigma_f*6)
    
    alphas = np.linspace(0,1,1000)
    for R in [-3+sing, sing, 3+sing]:
        poly = poly_for_alpha(alphas*R,R)
        plt.plot(alphas,poly,label=R)


def polynomial(R):
    
    A = (3 * f_p / (4*np.pi))**(2/3) * np.pi
    B = mu * b /(sigma_f*R)

    # print('B',B)

    poly = np.poly1d([6-B,-12+6*B,8-6/A-12*B,8*B+12/A,-8/A])
    return poly



R_array = np.linspace(19.5,122,2000)
def root_behaviour(R_array):
    sing = mu*b/(sigma_f*6)
    alpha_list = []
    for R in R_array:

        poly = polynomial(R)
        
        for x in poly.r:
            if 0 < x.real < 1:
                alpha_list.append(x.real)
    
    return np.array(alpha_list)
alpha_array = root_behaviour(R_array)

# plt.plot(np.linspace(0.01,100,1000),(alpha_array/(2-alpha_array))**3)

# plt.plot(np.linspace(0.01,100,1000),alpha_array**3 * (3*alpha_array-2))


def asympt_He_M_func_alpha1(R_array):
    A = 18.573
    B = -7.100
    C = 5.375

    p = (2 * gamma + 0.5 * mu * b) / R_array

    

    HeM = f_p / (n_mh * A * p**(-1))
    HeM = f_p / (n_mh * func_v_a(p))
    return HeM

def asympt_He_M_func_alpha0(R_array):
    A = 18.573
    B = -7.100
    C = 5.375
    new_A = (3 * f_p / (4*np.pi))**(2/3) * np.pi
    r = sigma_f / (new_A * mu * b) * R_array**2
    

    p = (2 * gamma + mu * b) / r

    v_a = A * p**(-1/3)

    HeM = f_p / (n_mh * v_a) * (r/R_array)**3
    return HeM

HeM = He_M_func(alpha_array*R_array,R_array)

R_array1 = np.linspace(50,150,1000)
# HeMasymp = asympt_He_M_func_alpha1(R_array1)
# a_array1 = root_behaviour(R_array1)
# plt.plot(HeMasymp,a_array1,label=r'asymp $\alpha\approx 1$',linestyle='dashed',color='blue')

R_array2 = np.linspace(0.1,20,1000)

HeMasymp2 =asympt_He_M_func_alpha0(R_array2)
 
# plt.plot(HeM,alpha_array,color='blue')
# plt.plot(R_array,func_p(alpha_array*R_array,R_array))


a_array2 = root_behaviour(R_array2)

# plt.plot(HeMasymp2,a_array2,label=r'asymp $\alpha\ll 1$',linestyle='dashed',color='blue')
# plt.xlabel('HeM')
# plt.ylabel(r'$\alpha$')
# plt.plot(x,poly(x),label='poly')

# plot_polynomial()





r_solution = alpha_array * R_array
He_M = HeM

He_interpolate = sp.interpolate.CubicSpline(r_solution, He_M)
R_function_r = sp.interpolate.CubicSpline(r_solution,R_array)
R_fracture = R_function_r(r_solution)
He_crit = He_interpolate(r_solution)

# print(R_fracture)


plt.plot(HeM,R_array)
plt.savefig('Helium/fracture_curve.eps',format='eps',bbox_inches='tight')

# plt.show()


v_total = 1000
He_array = np.linspace(0.01,0.485,v_total)

V_f = np.zeros(v_total)
V_t = np.zeros(v_total)
dx = R_array[1] - R_array[0]
for i in range(v_total):
    
    for j in range(0,len(R_array)):
        if R_array[j] > R0: 
            V_t[i] += dx * R_array[j]**3 * R_dist_func(R_array[j])    
            if He_array[i] > He_M[j]:
                V_f[i] += dx * R_array[j]**3 * R_dist_func(R_array[j])
V_f = V_f / V_t
# print(V_f)


Rp = 0.5 * 1e1 * 1e3  # 1 um
L = 1000
D_cl = 2 * R0 / ( 1 - (V_f/V_c)**(0.675) )
# D_cl = 16.67 / (V_c - V_f)**0.675
RF = V_f * (1- (1-D_cl/Rp)**3)
# RF = V_f * (D_cl / L)


print('V_f',V_f)

# print(RF)

for i in range(len(RF)):
    if RF[i] < 0:
        RF[i] = 0.2
    elif RF[i] > 0.2:
        RF[i] = 0.2


plt.cla()
plt.plot(He_array,RF)

np.savetxt('He_array.txt',He_array)
np.savetxt('RF.txt',RF)
print('saved')

# plt.savefig('Helium/percolation_curve.eps',format='eps',bbox_inches='tight')

# plt.legend()
plt.show()


# At some value of He_M... Find the R solution at the two points. Compute sum over f(R) dR 





# plt.plot(R_array,(r_solution/R_array)**3)
# plt.plot(R_array,1/v_a)

# A = (3 * f_p / (4*np.pi))**(2/3) * np.pi

# plt.plot(R_array,r_solution/R_array)
# plt.plot(R_array,((1/A)**0.5 * np.ones(len(R_array))),label=r'$\alpha=1/\sqrt{A}$')

# plt.plot(R_array,r_solution)
# plt.plot(R_array,sigma_f / (A*mu*b) * R_array**2,label='approx small R')

# print(1/A,'1/A')

# plt.plot(R_array,1/np.sqrt(np.pi) * (4*np.pi/(3*f_p))**(1/3) * np.ones(len(R_array)))

# plt.xlabel(r'Bubble spacing ')
# plt.ylabel(r'$\alpha = $ Bubble radius / Bubble spacing')

# plt.savefig('Helium/limitRzero.eps',format='eps',bbox_inches='tight')


# plt.plot(r_solution,v_a)

# plt.plot(R_array,((f_p * 3 / (4*np.pi))**(-2/3) / (np.pi))**0.5 * R_array,label='asympt')


# plt.plot(He_M,r_solution/5)
# plt.xlabel('He / M')
# plt.ylabel('bubble diameter (nm)')
# plt.xlim(0.35,0.6)
# plt.ylim(0,40)

# plt.legend()
# plt.show()

# plt.plot(He_array,RF)
# plt.xlabel('He / M')
# plt.ylabel('RF')

# plt.show()
