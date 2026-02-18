import numpy as np
from matplotlib import pyplot as plt


data = np.array([0.015338994615134327, 25.543525040150314,
0.017196930731522994, 20.622828471558883,
0.07602527632863428, 9.76660333588756,
0.09331772806953088, 8.57815403025182,
0.17286888429361685, 5.143176545918315,
0.24919541919027577, 4.214628360292654,
0.2951588692831726, 3.5663451142578246,
0.35028813755025345, 2.642625463171928,
0.5014002750165325, 2.425971217735416])

data_formatted = data.reshape((2,9),order='F')

inverse_square_root_p = data_formatted[0]

inversed = 1/inverse_square_root_p

squared = inversed**2

plottable = squared *0.00133322

plt.scatter(plottable**0.3333,data_formatted[1])

import scipy as sp

print(sp.stats.linregress(np.log(plottable),np.log(1/data_formatted[1])))



plt.show()