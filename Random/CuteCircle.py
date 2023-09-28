# Using M1 and M2 complex matrices 500x500
# when calculating eigenvalues of M1(1-t)*M2*t, the eigenvalues are complex
# for times inside the interval [0,1]

import numpy as np
import matplotlib.pyplot as plt

msize = 150
ts = 1000

def Mod_GUE(M1, M2, t):
    # We generate a random matrix
    M = (1-t)*M1 + t*M2
    # We compute the eigenvalues
    eig = np.linalg.eigvals(M)
    # return unique eigenvalues
    return eig

# for t in np.linspace(0,1,0.01):
times = np.linspace(0,1,ts)
# list of complex values
eig = [[complex(0,0) for i in range(msize)] for j in range(ts)]

# create 2 random matrices
M1 = np.random.randn(msize,msize) + 1j*np.random.randn(msize,msize)
M2 = np.random.randn(msize,msize) + 1j*np.random.randn(msize,msize)

for i in range(ts):
    print(i)
    eig[i] = Mod_GUE(M1, M2, times[i])

# create len(times) colors in a color shade from red to blue
colors = plt.cm.RdBu(np.linspace(0,1,ts))

# plot the eigenvalues for each time with colorcode
for i in range(ts):
    print(i)
    plt.scatter(eig[i].real, eig[i].imag, color=colors[i], s=0.05)
    # for x,y in zip(eig[i].real, eig[i].imag):
    #     plt.scatter(x,y, color=colors[i], s=0.1)

# remove axis
plt.axis('off')

# save the plot with high resolution
plt.savefig('CuteCircleWhite.png')

# black background
plt.style.use('dark_background')

# save the plot
plt.savefig('CuteCircleBlack.png')

