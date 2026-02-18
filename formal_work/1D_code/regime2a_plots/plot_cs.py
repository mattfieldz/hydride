from matplotlib import pyplot as plt
import numpy as np
import scipy as sp


Nstar = 8.01e-2
R=8.31

T = np.linspace(273.15+75,273.15+200,2)
P = 0.1


def Csstar(T):
    return 2.5e-1 * np.exp(-52000/(R*T))

def Drefstar(T):
    return 1.9e-2 * np.exp(-48387/(R*T))  
def kstar(T):
    return 10.4/Nstar * np.exp(1592/T)
def Jstar(T,P):
    return 1597 * np.sqrt(P) * np.exp(-88000/(R*T))

Lrefstar = np.sqrt( Drefstar(T) / (Nstar * kstar(T)))

Crefstar = Jstar(T,P) * Lrefstar / Drefstar(T)
print(T-273.15)
print(Csstar(T)/Crefstar)