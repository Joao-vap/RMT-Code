# load Memory.txt

import numpy as np
import matplotlib.pyplot as plt

Memory = np.loadtxt('Memory_full_med.txt')

plt.plot(Memory)
plt.show()

print(Memory.shape)

memory_max = Memory[:,-1]
memory_min = Memory[:,0]
memory_med_max = Memory[:,int(Memory.shape[1]/2)-1]
memory_med_min = Memory[:,int(Memory.shape[1]/2)]

# plot all lines in one plot
fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)

axs.set_xlabel('Time')
axs.set_ylabel('Lambdas')

axs.plot(memory_max, label='max')
axs.plot(memory_min, label='min')
axs.plot(memory_med_max, label='med min')
axs.plot(memory_med_min, label='med max')

axs.legend()

plt.show()
