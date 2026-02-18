import numpy as np

pure = [19,16,]
applied = [10,20,]
applied_inflated = [pure[i]+1 if pure[i] > applied[i] else applied[i] for i in range(len(pure))]

print(applied_inflated)

