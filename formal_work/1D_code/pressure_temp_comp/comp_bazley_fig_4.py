import numpy as np
from matplotlib import pyplot as plt
plt.style.use('formal_work/minorticks.mplstyle')

data_75C = np.array([0.060519175702996524, 0.00018957356524063754,
0.060519175702996524, 0.0000807992644838296,
0.1000000000000001, 0.0001291549665014884,
0.09905690937711491, 0.0001704046656425085,
0.1000000000000001, 0.0002020949938191079,
0.29733841963526814, 0.0006124026157621058,
0.29733841963526814, 0.00037503696379226933,
0.49599037752380687, 0.0012376350284724815,
0.49599037752380687, 0.0008613609622059082,
0.9905690937711485, 0.004447829676127635,
0.9905690937711485, 0.0030302710828663993,

])

data_100C = np.array([0.10095206950107322, 0.000774263682681127,
0.1000000000000001, 0.0005868383921470198,
0.5054796821191241, 0.004447829676127635,
1, 0.008431909292866269,
1.0095206950107323, 0.010660504989847927,
])

data_130C = np.array([0.060519175702996524, 0.0009380418666398144,
0.05994842503189415, 0.0015317404637020815,
0.1000000000000001, 0.0015647481416580207,
0.1000000000000001, 0.002296736176338636,
0.3001692880435886, 0.014994290682819923,
0.3001692880435886, 0.01376857164852759,
0.5007125506364689, 0.023462288481422636,
0.5007125506364689, 0.022482875308090036,
1.0095206950107323, 0.03162277660168381,
1, 0.016329126856448006,
])




# These are p values and corresponding 1/t

data_formatted_75C = data_75C.reshape((2,11),order='F')
data_formatted_100C = data_100C.reshape((2,5),order='F')
data_formatted_130C = data_130C.reshape((2,10),order='F')

P=0.1
T = 273.15 + 75 

R=8.31
Nstar = 8.01e-2

alpha_b = 0.8
alpha_ck = 0.99


def Drefstar(T):
    return 1.9e-2 * np.exp(-48387/(R*T))  
def kstar(T):
    return 10.4/Nstar * np.exp(1592/T)
def Jstar(T,P):
    return 1597 * np.sqrt(P) * np.exp(-88000/(R*T))
def alpha_c(P):
    return alpha_ck - 0.1 * P



def tstar_b(T,P):
    
    return np.log(1/alpha_c(P)) / (2 * Jstar(T,P)) * np.sqrt(3*Drefstar(T)*Nstar/kstar(T)) * (1 + (np.log(alpha_b)/np.log(alpha_c(P)))**2)


p_arr = np.linspace(0.1,1,4)




def Castarck(T,P):
    return Nstar * 4.13e-6 * np.exp(-894/T) * np.sqrt(1e5 * P)

def condon_kirkpatrick_displacement(T,P,time_seconds):
    return np.sqrt(Drefstar(T)*kstar(T)*(Castarck(T,P))**(2)/(3*Nstar))/(1-alpha_ck)*time_seconds

displacement = 10 * 1e3 * 1e-7



def tstar_ck(T,P):
    return displacement*np.sqrt(3*Nstar/(Drefstar(T) * kstar(T))) * (1-alpha_ck) / Castarck(T,P)

print(tstar_ck(75+273.15,0.1))
print(tstar_b(75+273.15,0.1))

def ttotal(T,P):
    return tstar_ck(T,P) + tstar_b(T,P)

plt.scatter(data_formatted_75C[0,:],data_formatted_75C[1,:],marker='x',color='red',label='75C')
plt.scatter(data_formatted_100C[0,:],data_formatted_100C[1,:],marker='x',color='blue',label='100C')
plt.scatter(data_formatted_130C[0,:],data_formatted_130C[1,:],marker='x',color='green',label='130C')

a_75, b_75 = np.polyfit(np.log(data_formatted_75C[0,:]),np.log(data_formatted_75C[1,:]), deg=1)
a_100, b_100 = np.polyfit(np.log(data_formatted_100C[0,:]),np.log(data_formatted_100C[1,:]), deg=1)
a_130, b_130 = np.polyfit(np.log(data_formatted_130C[0,:]),np.log(data_formatted_130C[1,:]), deg=1)
regress_p = np.linspace(0.01,10,10)


plt.plot(regress_p,np.exp(b_75)*regress_p**a_75,color='red')
plt.plot(regress_p,np.exp(b_100)*regress_p**a_100,color='blue')
plt.plot(regress_p,np.exp(b_130)*regress_p**a_130,color='green')

plt.scatter(p_arr,1/ttotal(75+273.15,p_arr),color='red')
plt.scatter(p_arr,1/ttotal(100+273.15,p_arr),color='blue')
plt.scatter(p_arr,1/ttotal(130+273.15,p_arr),color='green')

# plt.scatter(p_arr,1/tstar_b(75+273.15,p_arr),color='red',marker='^')
# plt.scatter(p_arr,1/tstar_b(100+273.15,p_arr),color='blue',marker='^')
# plt.scatter(p_arr,1/tstar_b(130+273.15,p_arr),color='green',marker='^')

plt.plot(p_arr,1/tstar_b(75+273.15,p_arr),color='red',linestyle='dashed')
plt.plot(p_arr,1/tstar_b(100+273.15,p_arr),color='blue',linestyle='dashed')
plt.plot(p_arr,1/tstar_b(130+273.15,p_arr),color='green',linestyle='dashed')



plt.xscale("log")
plt.yscale("log")

plt.xlabel('pressure (bar)')
plt.ylabel('1 / initation time (seconds)')

plt.xlim(0.1*min(p_arr),10*max(p_arr))
plt.ylim(1e-5,1e-1)

plt.legend()
plt.show()