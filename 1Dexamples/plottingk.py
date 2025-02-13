import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import linregress
vec = []


m_vec = []
values = [1e4,1e5,5e5,1e6,5e6,1e7]
j = 0
for i in values:
    
    vec.append(np.loadtxt(f"1Dexamples/data/k{i}.dat"))
    m, b = np.polyfit(np.linspace(0,10,1001), vec[j], 1)
    m_vec.append(m)
    j += 1
print(m_vec)


# m1,b1 = np.polyfit(values,m_vec,1)

# a = linregress(values,m_vec)
# print(a)

m_array = np.array(m_vec)
k_array = np.array(values)

m_log = np.log(m_array)
k_log = np.log(k_array)

# plt.plot(values,m_vec)

print(np.polyfit(k_log,m_log,1))

plt.plot(k_log,m_log)

plt.show() 