import numpy

t = 0
dt = 1/256
for i in range(1,1000):
    t += dt
    # t = i*dt
print(t)