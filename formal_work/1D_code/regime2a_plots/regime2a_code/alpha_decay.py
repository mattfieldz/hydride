import numpy as np
import scipy
alpha_c = 0.98

tau_c = -np.sqrt(3) * np.log(alpha_c)

print(tau_c)

alpha_initial = alpha_c
gamma_initial = 0
t_initial = 0


def f(x):
    return 1/(np.sqrt(3)*(1+np.sqrt(3)*x))

tmax = 1
dt = 0.01



alpha = alpha_initial
gamma = gamma_initial
t = t_initial

t_list = [0]
alpha_list = [alpha_c]
gamma_list = [gamma]

while t < tmax:
    t += dt


    gamma_new = gamma + dt * (np.exp(np.sqrt(3)*gamma)/(tau_c * (1+np.sqrt(3)*gamma)))

    alpha_new = alpha - dt * f(gamma) * alpha
    # gamma_new = gamma - 1/np.sqrt(3) * np.log(np.sqrt(3)*(1+np.sqrt(3) * gamma) * (alpha - alpha_c))


    t_list.append(t)
    alpha_list.append(alpha_new)
    gamma_list.append(gamma_new)
    alpha = alpha_new
    gamma = gamma_new


# print(alpha_list)
print(gamma_list)