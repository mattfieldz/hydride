import numpy as np
from matplotlib import pyplot as plt
from copy import copy

plt.style.use('formal_work/minorticks.mplstyle')
# plt.margins(0)


# exp_1 = np.array([628.0405405405407, 74.14187643020614,
# 1062.837837837838, 183.29519450800922,
# 1493.9189189189187, 362.47139588100686,
# 2058.7837837837837, 630.2059496567506,
# 2460.1351351351354, 813.5011441647596,
# ])

# exp_2 = np.array([628.0405405405407, 61.784897025171745,
# 1062.837837837838, 158.58123569794066,
# 1493.9189189189187, 331.57894736842104,
# 2058.7837837837837, 591.0755148741418,
# 2460.1351351351354, 768.1922196796338,
# ])

# exp_3 = np.array([583.4459459459462, 24.71395881006879,
# 1018.2432432432435, 55.606407322654604,
# 1449.3243243243242, 109.1533180778033,
# 2017.9054054054054, 191.53318077803215,
# 2415.5405405405404, 294.5080091533181,
# ])

# exp_4 = np.array([583.4459459459462, 28.83295194508014,
# 1018.2432432432435, 74.14187643020614,
# 1453.0405405405404, 146.22425629290626,
# 2017.9054054054054, 325.4004576659039,
# 2415.5405405405404, 523.1121281464532,
# ])

# exp_5 =np.array([966.2162162162161, 12.356979405034394,
# 1367.5675675675675, 16.475972540045746,
# ])

# new data ... 

exp_1 = np.array([631.9702602230484, 75.32956685499084,
1055.7620817843865, 188.32391713747666,
1494.4237918215615, 362.5235404896425,
2059.4795539033457, 630.8851224105463,
2460.96654275093, 814.5009416195858,
3947.955390334573, 1445.3860640301323,
5732.342007434944, 2109.227871939737,
])

exp_2 = np.array([631.9702602230484, 61.20527306968006,
1063.197026022305, 160.0753295668551,
1494.4237918215615, 329.56685499058403,
2059.4795539033457, 588.5122410546142,
2468.401486988848, 762.71186440678,
3955.3903345724902, 1412.4293785310738,
5739.776951672862, 2066.8549905838045,
])

exp_3 = np.array([1018.5873605947957, 61.20527306968006,
1457.2490706319702, 108.28625235404934,
2022.3048327137546, 188.32391713747666,
2423.791821561338, 291.9020715630886,
3910.7806691449814, 1012.2410546139363,
5695.167286245354, 1699.6233521657255,
])

exp_4 = np.array([587.3605947955389, 28.248587570621567,
1018.5873605947957, 70.62146892655392,
1457.2490706319702, 141.24293785310783,
2022.3048327137546, 324.8587570621471,
2423.791821561338, 522.5988700564976,
3910.7806691449814, 1228.8135593220343,
5695.167286245354, 1864.4067796610173,
])

exp_5 = np.array([973.9776951672864, 14.124293785311238,
1375.4646840148696, 18.832391713748166,
2862.453531598513, 42.37288135593235,
4646.840148698884, 84.7457627118647,
])




e1 = exp_1.reshape((2,7),order='F')
e2 = exp_2.reshape((2,7),order='F')
e3 = exp_3.reshape((2,6),order='F')
e4 = exp_4.reshape((2,7),order='F')
e5 = exp_5.reshape((2,4),order='F')


experiments = [e1,e2,e3,e4,e5]

original_experiments_data = [np.copy(e1),np.copy(e2),np.copy(e3),np.copy(e4),np.copy(e5)]



# e1 87.9% tritium 19.5% load 
# e2 87.9% tritium 19.8% load
# e3 39.0% tritium 59.5% load
# e4 42.5% tritium 55.6% load
# e5 15.2% tritium 54.7% load

ratios = np.array([0.879,0.879,0.39,0.425,0.152])


# pV = nRT 
# p = 100,000 (1atm)
# V in cm^3
# R = 8.31
# T=273.15

v_factor_to_convert_to_mol_of_He = 101325 / (8.31*273.15*100**3)


for e in experiments:
    e[1,:] = e[1,:] * v_factor_to_convert_to_mol_of_He
    
# This has converted the cm^3 into moles of helium released

# molar mass of tritium is 3.016 g / mol

# moles of tritium is grams / molar mass (varies experiment to experiment)

grams_of_tritium = np.array([0.703,0.696,0.915,0.938,0.336])

# not diatomic
moles_of_tritium = (grams_of_tritium / 3.016 ) 

# decay constant = log(2) / half life
# half life of tritium is 12.33 years

decay_constant = np.log(2) / (12.32*365)


index = 0
for e in experiments:
    
    e[0,:] = moles_of_tritium[index] * (1-np.exp(-decay_constant * e[0,:]))
    index+=1
# This has converted into moles of Helium generated

for e in experiments:
    e[1,:] = (e[1,:] / e[0,:]) 
# This has converted to release fraction



U_eff = np.copy(moles_of_tritium)
U_eff = (1 / ratios * (U_eff)) / 3



index=0
for i in range(5):
    index+=1

    # Northall style
    # plt.scatter((experiments[i][0,:]* (1-experiments[i][1,:]))/U_eff[i],experiments[i][1,:],label=index)

    # Just the helium generated on the x axis
    # plt.scatter((experiments[i][0,:])/U_eff[i],experiments[i][1,:],label=index,marker='x')
    
    print((experiments[i][0,:])/U_eff[i],experiments[i][1,:])

    plt.scatter(original_experiments_data[i][0,:],original_experiments_data[i][1,:],marker='x')

    np.savetxt(f'Helium/diffusion_model/data/original_experiments_{i}.txt',original_experiments_data[i])
    # plt.scatter((original_experiments_data[i][0,:]),experiments[i][0,:]* (1-experiments[i][1,:])/U_eff[i])
    # print(original_experiments_data[i][0,:])
    # plt.scatter((experiments[i][0,:]/(U_eff[i])),experiments[i][1,:])

t_vals = np.linspace(0,6000,100)
He_production =  (1 - np.exp(-decay_constant * t_vals))

He_production = He_production / v_factor_to_convert_to_mol_of_He / 3.016 * 0.7 

plt.plot(t_vals,He_production)




# plt.savefig('Helium/northall_rescaled.eps',format='eps',bbox_inches='tight')

He_array = np.loadtxt('He_array.txt')
RF = np.loadtxt('RF.txt')

print(len(He_array))
# plt.plot(He_array[0:600],RF[0:600])

# plt.savefig('Helium/northall_with_model.eps',format='eps',bbox_inches='tight')


plt.xlabel(r'Time (days)')
plt.ylabel(r'He released $(cm^3)$')
plt.legend()
plt.show()