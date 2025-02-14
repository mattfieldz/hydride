import numpy as np
import scipy as sp
from matplotlib import pyplot as plt


total_time = 1
a = 0.75
sigma = 0.6


dx_list = [0.1,0.05,0.02,0.01,0.005,0.002,0.001]
error_list = []
for dx in dx_list:
    t = 0 

    n = int(2/dx)

    dt = dx * (sigma)/a

    x = np.linspace(-1,1,n)

    u_prev = np.sin(2*np.pi*x)
    u_current = np.sin(2*np.pi*(x-a*dt))
    u_update = np.zeros(n)
    # fig,axis = plt.subplots(2,2)

    error = []
    


    while t < total_time:

        

# leap frog
        u_update[0] = u_prev[0] - sigma * (u_current[1]-u_current[n-1])
        
        for i in range(1,n-1):
            u_update[i] = u_prev[i] - sigma * (u_current[i+1]-u_current[i-1])
        
        u_update[n-1] = u_prev[n-1] - sigma * (u_current[0]-u_current[n-2])
        

# angle derivative
        # u_update[0] = u_prev[n-1] + (1-2*sigma) * (u_current[0]-u_current[n-1])

        # for i in range(1,n):
        #     u_update[i] = u_prev[i-1] + (1-2*sigma) * (u_current[i]-u_current[i-1])



        t += dt

        u_prev = np.copy(u_current)
        u_current = np.copy(u_update)



        analytic = np.sin(2*np.pi*(x-a*t))

        # axis[0,0].cla()
        # axis[0,0].plot(x,u_current,label='approx')
        # axis[0,0].plot(x,analytic,label='analytic')
        
        error.append(np.linalg.norm(analytic-u_current))
    
        # leg = axis[0,0].legend()

        # fig.suptitle(t)
        # plt.pause(0.01)
        
    error_list.append(np.max(error))
    print(dx)
    
print(dx_list,error_list)

m,b = np.polyfit(np.log(dx_list),np.log(error_list),1)

print('m = ',m,'b = ',b)

plt.plot(np.log(dx_list),np.log(error_list))
plt.show()