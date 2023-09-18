# open file Memory.txt and plot histogram of first column

import numpy as np
import matplotlib.pyplot as plt
from TracyWidom import TracyWidom

Memory = np.loadtxt('Memory_t_cupula.txt')

t = 0.9840
n_particles = 1000

def sigma_t(t):
    return 2*t*(1-t)

aux_memory = [Memory[:,499], Memory[:,500]]

# subtract median
aux_memory[0] = aux_memory[0] - np.mean(aux_memory[0])
aux_memory[1] = aux_memory[1] - np.mean(aux_memory[1])

# plot histogram for (x-med)/npart**(-2/3) with subplot for each column
fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

data1 = n_particles**(3/4)*(aux_memory[0])/(sigma_t(t)*np.sqrt(n_particles))
data2 = n_particles**(3/4)*(aux_memory[1])/(sigma_t(t)*np.sqrt(n_particles))

axs[0].set_title('Lower')
axs[1].set_title('Upper')

# set bins
bins1 = np.linspace(min(data1), max(data1), 100)
bins2 = np.linspace(min(data2), max(data2), 100)

# get histogram normalized to 1
hist1, bins1 = np.histogram(data1, bins=bins1, range=(-60,60), density=True)
hist2, bins2 = np.histogram(data2, bins=bins2, range=(-60,60), density=True)

# get bin centers
bin_centers1 = 0.5*(bins1[1:] + bins1[:-1])
bin_centers2 = 0.5*(bins2[1:] + bins2[:-1])

# normalize histogram to 1
hist1 = hist1/np.sum(hist1)
hist2 = hist2/np.sum(hist2)

print(np.sum(hist1))

# plot histogram in axes
axs[1].plot(bin_centers1, hist1)
axs[0].plot(bin_centers2, hist2)

# plot median
m1 = np.median(data1)
axs[1].axvline(m1, color='k', linestyle='dashed', linewidth=1)
m2 = np.median(data2)
axs[0].axvline(m2, color='k', linestyle='dashed', linewidth=1)

# plot tracy-widom distribution centered at 0
x = np.linspace(-4,4,100)
tw = TracyWidom(2)
pdf = tw.pdf(x)

# center
x = x + 1.7710868074116265

# normalize pdf to 1
pdf = pdf/np.sum(pdf)

axs[0].plot(x, pdf)
axs[1].plot(-x, pdf)

# set labels
axs[0].set_xlabel(r'$N^{\frac{3}{4}} \left( \frac{\lambda_{min}}{\sigma \sqrt{N}}  \right)$')
axs[1].set_xlabel(r'$N^{\frac{3}{4}} \left( \frac{\lambda_{max}}{\sigma \sqrt{N}}  \right)$')
axs[0].set_ylabel('Probability')

plt.show()