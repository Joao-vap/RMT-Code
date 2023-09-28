import numpy as np
import matplotlib.pyplot as plt

Memory = np.loadtxt('Memory_t_cupula.txt')

plt.hist(Memory[:,499], bins=30, range=(-0.04,0.04), density=True)
plt.hist(Memory[:,500], bins=30, range=(-0.04,0.04), density=True)

plt.show()