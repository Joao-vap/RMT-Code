import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#plot two dimensional data

files = ['../Complex/testa02t04.dat']

data = []

beta = 2
N = 100

for file in files:
    dataaux = []
    with open(file) as f:
        for line in f:
            aux = line.split(" ")
            # filter empty string and \n
            aux = list(filter(lambda x: x != '\n' and x != '', aux))
            # make a list of complex numbers
            aux = np.array([float(i) for i in aux])
            # change shape to (N, 2)
            aux = aux.reshape(-1, N)
            for i in range(aux.shape[1]):
                dataaux.append(complex(aux[0, i], aux[1, i]))
    data.append(dataaux)

# plot two dimensional data
fig, ax = plt.subplots(1, 1, tight_layout=True)
# set proportion
plt.gca().set_aspect('equal')

# scatter
# ax.scatter([i.real for i in data[0]], [i.imag for i in data[0]], color='grey', s=1)

# #heatmap 2d greyscale

data = data[0]

# filter any complex number with abs > 1
data = list(filter(lambda x: abs(x) < 1, data))

#grid
space = 0.01
minx = min([i.real for i in data]) - space
maxx = max([i.real for i in data]) + space
miny = min([i.imag for i in data]) - space
maxy = max([i.imag for i in data]) + space
s = 500

Z = np.zeros((s, s))
for k in range(len(data)):
    part = data[k]
    coordx = int((part.real - minx) / (maxx - minx) * s)
    coordy = int((part.imag - miny) / (maxy - miny) * s)
    Z[coordy, coordx] += 1

#normalize
Z = Z / np.sum(Z)
ax.imshow(1 - Z, cmap='gray', interpolation='nearest')

ax.axis('off')

# ax[0,1].set_ylabel('a = 5')
# ax[1,1].set_ylabel('a = 3')
# ax[2,1].set_ylabel('a = 2')
# ax[3,1].set_ylabel('a = 1')

# ax[]

plt.savefig('tcomplexplot.png')

plt.show()