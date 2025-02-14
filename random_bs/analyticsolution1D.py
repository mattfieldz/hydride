import numpy as np
from matplotlib import pyplot as plt

length = 10
time_length = 3

n = 200

u = np.zeros(n)

t = 2

for i in range(n):
    c = 0
    for j in range(50):
        c = c + 1/j * np.exp(-3*j**2*np.pi**2*t) * np.sin(k*np.pi*x)

    u[i] = 1 - 1/np.pi * 2