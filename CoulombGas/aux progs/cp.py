import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#plot two dimensional data

# files = ['../Complex/testa02t001.dat', '../Complex/testa02t003.dat', '../Complex/testa02t005.dat', 
#          '../Complex/testa01t001.dat', '../Complex/testa01t003.dat', '../Complex/testa01t005.dat', '../Complex/testa01t007.dat',
#          '../Complex/testa01t009.dat', '../Complex/testa01t011.dat',  '../Complex/testa01t013.dat', 
#          '../Complex/testa00t001.dat', '../Complex/testa00t003.dat', '../Complex/testa00t005.dat', '../Complex/testa00t007.dat',
#          '../Complex/testa00t009.dat', '../Complex/testa00t011.dat', '../Complex/testa00t013.dat', '../Complex/testa00t015.dat',
#          '../Complex/testa00t017.dat', '../Complex/a00t019.dat', '../Complex/testa00t021.dat', '../Complex/testa00t023.dat',
#          '../Complex/testa00t025.dat', '../Complex/a00t027.dat',
#          '../Complex/testam01t001.dat', '../Complex/testam01t003.dat', '../Complex/testam01t005.dat', '../Complex/testam01t007.dat',
#          '../Complex/testam01t009.dat', '../Complex/testam01t011.dat', '../Complex/testam01t013.dat', '../Complex/testam01t015.dat',
#          '../Complex/testam01t017.dat', '../Complex/testam01t019.dat', '../Complex/testam01t021.dat', '../Complex/testam01t023.dat',
#          '../Complex/testam02t001.dat', '../Complex/testam02t003.dat', '../Complex/testam02t005.dat', '../Complex/testam02t007.dat',
#          '../Complex/testam02t009.dat', '../Complex/testam02t011.dat', '../Complex/testam02t013.dat', '../Complex/testam02t015.dat',
#          '../Complex/testam03t001.dat', '../Complex/testam03t003.dat', '../Complex/testam03t005.dat', '../Complex/testam03t007.dat',
#          '../Complex/testam03t009.dat', '../Complex/testam03t011.dat',
#          '../Complex/testam04t001.dat', '../Complex/testam04t003.dat', '../Complex/testam04t005.dat', '../Complex/testam04t007.dat',
#          '../Complex/testam05t001.dat']


files = ['../Complex/testa02t001.dat',  '../Complex/testa02t005.dat', 
         '../Complex/testa01t001.dat',  '../Complex/testa01t005.dat',  '../Complex/testa01t009.dat',  '../Complex/testa01t013.dat', 
         '../Complex/testa00t001.dat',  '../Complex/testa00t005.dat',  '../Complex/testa00t009.dat',  '../Complex/testa00t013.dat', '../Complex/testa00t017.dat','../Complex/testa00t021.dat', '../Complex/testa00t025.dat',
         '../Complex/testam01t001.dat', '../Complex/testam01t005.dat', '../Complex/testam01t009.dat', '../Complex/testam01t013.dat','../Complex/testam01t017.dat','../Complex/testam01t021.dat',
         '../Complex/testam02t001.dat', '../Complex/testam02t005.dat', '../Complex/testam02t009.dat', '../Complex/testam02t013.dat',
         '../Complex/testam03t001.dat', '../Complex/testam03t005.dat', '../Complex/testam03t009.dat',
         '../Complex/testam04t001.dat', '../Complex/testam04t005.dat',
         '../Complex/testam05t001.dat']

fig, ax = plt.subplots(8, 7, figsize=(20, 6), sharex=True, sharey=True, layout='tight')

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
            if len(aux) == 6*N:
                aux = aux.reshape(-1, 3*N)
            else:
                aux = aux.reshape(-1, N)
            for i in range(aux.shape[1]):
                dataaux.append(complex(aux[0, i], aux[1, i]))
    data.append(dataaux)

# reshape data to (...)
sizes = [2, 4, 7, 6, 4, 3, 2, 1]
dataaux = [[] for i in range(len(sizes))]
for d in range(len(dataaux)):
    dataaux[d] = data[sum(sizes[:d]):sum(sizes[:d+1])]

minx = -1.4
maxx = 1.4
miny = -1.4
maxy = 1.4

for j in range(len(dataaux)):
    for i in range(len(dataaux[j])):

        print(i, j)
        
        dat = dataaux[j][i]

        # filter any complex number with abs > 1
        dat = list(filter(lambda x: abs(x) < 1.4, dat))

        #grid
        #space = 0.1
        # minx = min([i.real for i in dat]) - space
        # maxx = max([i.real for i in dat]) + space
        # miny = min([i.imag for i in dat]) - space
        # maxy = max([i.imag for i in dat]) + space
        s = 1500

        Z = np.zeros((s, s))
        for k in range(len(dat)):
            part = dat[k]
            coordx = int((part.real - minx) / (maxx - minx) * s)
            coordy = int((part.imag - miny) / (maxy - miny) * s)
            Z[coordy, coordx] += 1

        # cut = 1 all that is bigger than 1
        Z = np.clip(Z, 0, 1)

        #Z = Z / np.sum(Z)
        ax[j][i].imshow(1-Z, cmap='gray', extent=[minx, maxx, miny, maxy], label='a = ' + str(j-4) + ', t = ' + str(i+1))
        #imshow(1 - Z, cmap='gray', interpolation='nearest')

# for i in range(8):
#     for j in range(14):
#         ax[i][j].axis('off')

# remoe tikc, labels and title
for i in range(8):
    for j in range(7):
        ax[i][j].set_xticks([])
        ax[i][j].set_yticks([])
        ax[i][j].set_xticklabels([])
        ax[i][j].set_yticklabels([])
        ax[i][j].set_title('')
        ax[i][j].set_xlabel('')
        ax[i][j].set_ylabel('')
        ax[i][j].spines['top'].set_visible(False)
        ax[i][j].spines['right'].set_visible(False)
        ax[i][j].spines['bottom'].set_visible(False)
        ax[i][j].spines['left'].set_visible(False)

# for bottom row set xlabel t_0  = 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27
for i in range(7):
    ax[7][i].set_xlabel('$t_0$ = ' + str(round((0.01 + 0.04*i)/2, 3)))

# for first line set ylabel a = 2, 1, 0, -1, -2, -3, -4, -5
for i in range(8):
    ax[i][0].set_ylabel('$t_1$ = ' + str(2-i))


# ax[0,1].set_ylabel('a = 5')
# ax[1,1].set_ylabel('a = 3')
# ax[2,1].set_ylabel('a = 2')
# ax[3,1].set_ylabel('a = 1')

# ax[]

plt.subplots_adjust(wspace=0, hspace=0)

plt.legend()

plt.savefig('tcomplexplot.png')

plt.show()