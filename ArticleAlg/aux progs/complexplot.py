import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#plot two dimensional data

files = ['../Complex/a5t0.dat', '../Complex/a5t05.dat', '../Complex/a5t1.dat', '../Complex/a5t2.dat',
         '../Complex/a3t0.dat', '../Complex/a3t05.dat', '../Complex/a3t1.dat', '../Complex/a3t2.dat', 
         '../Complex/a2t0.dat', '../Complex/a2t05.dat', '../Complex/a2t1.dat', '../Complex/a2t2.dat',]

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
fig, ax = plt.subplots(3, 4, tight_layout=True)
# set proportion
plt.gca().set_aspect('equal')

# scatter
# ax.scatter([i.real for i in data[0]], [i.imag for i in data[0]], color='grey', s=1)

# #heatmap 2d greyscale

#grid
l = 1.5
s = 200

for i in range(3):
    for t in range(4):
        print(i, t)
        adata = data[i * 4 + t]
        Z = np.zeros((s, s))
        for k in range(len(adata)):
            part = adata[k]
            coordx = int((part.real + l) / (2 * l) * s)
            coordy = int((part.imag + l) / (2 * l) * s)
            Z[coordx, coordy] += 1

        #normalize
        Z = Z / np.sum(Z)
        ax[i, t].imshow(1 - Z, cmap='gray', interpolation='nearest')

#turn off axis
for i in range(3):
    for j in range(4):
        ax[i, j].axis('off')

ax[0, 0].set_title('t = 0')
ax[0, 1].set_title('t = 0.5')
ax[0, 2].set_title('t = 1')
ax[0, 3].set_title('t = 2')

ax[0,0].set_ylabel('a = 5')
ax[1,0].set_ylabel('a = 3')
ax[2,0].set_ylabel('a = 2')

plt.savefig('complexplot.png')

plt.show()