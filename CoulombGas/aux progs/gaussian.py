#read file dataX.txt and plot histogram
import matplotlib.pyplot as plt
import numpy as np
from TracyWidom import TracyWidom
import math

files = ['../Gaussian/b1N10.dat', '../Gaussian/b1N100.dat',
         '../Gaussian/b2N10.dat', '../Gaussian/bN2N100.dat',
         '../Gaussian/b4N10.dat', '../Gaussian/b4N100.dat',]

data = []

for file in files:
    dataaux = []
    with open(file) as f:
        for line in f:
            aux = line.split(" ")
            # filter empty string and \n
            aux = list(filter(lambda x: x != '\n' and x != '', aux))
            # make a list of floats
            aux = list(map(float, aux))
            for i in aux:
                dataaux.append(i)
    data.append(dataaux)

# were doing two subplots in one figure
fig, axs = plt.subplots(3, 3, sharey=True, tight_layout=True)

# title for the whole figure
fig.suptitle('Ensembles Gaussianos - Densidade')

# for the first subplot
for i in range(3):

    if i == 0:
        beta = 1
    elif i == 1:
        beta = 2
    else:
        beta = 4

    for j in range(2):

        if j == 0:
            N = 10
            adata = [i / math.sqrt(2*beta) for i in data[2*i + j]]
        else:
            N = 100
            adata = [i / math.sqrt(2*beta) for i in data[2*i + j]]

        axs[i, j].hist(adata, bins=50, range=(-1.5, 1.5), density=True, color='gray', edgecolor='black',linewidth=0.5)
        axs[i, j].set_title("beta = " + str(beta) + ", N = " + str(N))
        axs[i, j].set_xlim(-1.5, 1.5)

files = ['../Gaussian/b1N100.dat',
         '../Gaussian/b2N100.dat',
         '../Gaussian/b4N100.dat',]

data = []

N = 100

for file in files:
    dataaux = []
    with open(file) as f:
        for line in f:
            aux = line.split(" ")
            # filter empty string and \n
            aux = list(filter(lambda x: x != '\n' and x != '', aux))
            # make a list of floats
            aux = list(map(float, aux))
            # for i in aux:
            #     dataaux.append(i)
            # append max value in the list
            dataaux.append(max(aux))
    data.append(dataaux)

for i in range(3):

    if i == 0:
        beta = 1
    elif i == 1:
        beta = 2
    else:
        beta = 4

#    adata = [(j - math.sqrt(2*beta)) * N**(4/6) for j in data[i]]
    adata = [(j * 2 / math.sqrt(2*beta)  - 2) * N**(2/3) for j in data[i]]

    # center the data
    #eliminate outliers
    adata = [i for i in adata if i < 5 and i > -5]
    adata = [i - np.mean(adata) for i in adata]

    axs[i, 2].hist(adata, bins=80, range=(-5, 5), color='gray', density=True, edgecolor='black',linewidth=0.5)
    axs[i, 2].set_title("Maior autovalor - beta = " + str(beta) + f", N = {N}")
    axs[i, 2].set_xlim(-4, 4)


axs[1, 0].set_ylabel("Frequência")
axs[2, 1].set_xlabel("Posição")

# plot the semi-circle for the axs[0,2], axs[1,2] and axs[2,2]
x = np.linspace(-3, 3, 1000)
y = 2/np.pi * np.sqrt(1 - x**2)

axs[0, 1].plot(x, y, color='red', alpha=0.8,label='Wigner Semi-Circle')
axs[1, 1].plot(x, y, color='red',alpha=0.8, label='Wigner Semi-Circle')
axs[2, 1].plot(x, y, color='red', alpha=0.8,label='Wigner Semi-Circle')

axs[0, 1].legend()
axs[1, 1].legend()
axs[2, 1].legend()

x = np.linspace(-4, 4, 1000)

# # We plot the theoretical density centered
y1 = TracyWidom(1).pdf(x - 1.2065335745820)
y2 = TracyWidom(2).pdf(x - 1.771086807411)
y4 = TracyWidom(4).pdf(x - 2.306884893241)
# y1 = TracyWidom(1).pdf(x)
# y2 = TracyWidom(2).pdf(x)
# y4 = TracyWidom(4).pdf(x)
axs[0, 2].plot(x, y1, label='Tracy-Widow', color='green', alpha=0.8)
axs[1, 2].plot(x, y2, label='Tracy-Widow', color='green', alpha=0.8)
axs[2, 2].plot(x, y4, label='Tracy-Widow', color='green', alpha=0.8)


# show labels
axs[0, 2].legend()
axs[1, 2].legend()
axs[2, 2].legend()

# We are going to compute the density of eigenvalues of a Random Matrix
# in the models of GOE, GUE and GSE

# number of eigenmvalues
N = 10

# We define the number of realizations
M = 1000

# Gaussian Orthogonal Ensemble
def GOE(N, escale = False):
    beta = 1
    # We generate a random matrix
    A = np.random.randn(N,N)
    # We make it symmetric
    A = (A + A.T)/2
    # We compute the eigenvalues
    eig = np.linalg.eigvalsh(A)
    # return unique eigenvalues
    if escale:
        return np.unique(eig)/(np.sqrt(N*beta))
    return np.unique(eig)

# Gaussian Unitary Ensemble
def GUE(N, escale = False):
    beta = 2
    # We generate a random matrix
    A = np.random.randn(N,N) + 1j*np.random.randn(N,N)
    # We make it Hermitian
    A = (A + A.T.conj())/2
    # We compute the eigenvalues
    eig = np.linalg.eigvalsh(A)
    # return unique eigenvalues
    if escale:
        return np.unique(eig)/(np.sqrt(N*beta))
    return np.unique(eig)

# Gaussian Symplectic Ensemble
def GSE(N, escale = False):
    beta = 4
    # We generate a random matrix with quaternions entries
    A = np.random.randn(N,N) + 1j*np.random.randn(N,N)
    B = np.random.randn(N,N) + 1j*np.random.randn(N,N)
    M = np.block([[A, -B.conj()], [B, A.conj()]])
    # We make it Hermitian
    M = (M + M.T.conj())/2
    # We compute the eigenvalues
    eig = np.linalg.eigvalsh(M)
    # return unique eigenvalues
    if escale:
        return np.unique(eig)/(np.sqrt(N*beta))
    return np.unique(eig)

# We compute the eigenvalues of the random matrices
eig_GOE = np.zeros((M,N))
eig_GUE = np.zeros((M,N))
eig_GSE = np.ones((M,2*N))*1000

for i in range(M):
    eig_GOE[i,:] = GOE(N, escale=True)
    eig_GUE[i,:] = GUE(N, escale=True)
    # complete GSE with 1000 as eigenvalues
    value = GSE(N, escale=True)
    eig_GSE[i,:len(value)] = value

# We compute the density of eigenvalues
x = np.linspace(-5,5,1000)
dx = x[1] - x[0]

density_GOE = np.zeros(x.shape)
density_GUE = np.zeros(x.shape)
density_GSE = np.zeros(x.shape)

eig_GOE = eig_GOE / np.sqrt(2)
eig_GUE = eig_GUE / np.sqrt(2)
eig_GSE = eig_GSE / np.sqrt(2)

for i in range(M):
    density_GOE[:-1] += np.histogram(eig_GOE[i,:], bins=x)[0]
    density_GUE[:-1] += np.histogram(eig_GUE[i,:], bins=x)[0]
    density_GSE[:-1] += np.histogram(eig_GSE[i,:], bins=x)[0]

density_GOE /= (M*N)/100
density_GUE /= (M*N)/100
density_GSE /= (M*N*2)/100

# # We plot the density of eigenvalues
# plt.plot(x, density_GOE, label='GOE', color='red', alpha=0.4)
# plt.plot(x, density_GUE, label='GUE', color='blue', alpha=0.4)
# plt.plot(x, density_GSE, label='GSE', color='green', alpha=0.4)

axs[0,0].plot(x, density_GOE, color='blue', alpha=0.4, label='Matrix Model')
axs[1,0].plot(x, density_GUE, color='blue', alpha=0.4, label='Matrix Model')
axs[2,0].plot(x, density_GSE, color='blue', alpha=0.4, label='Matrix Model')

#legend
axs[0,0].legend()
axs[1,0].legend()
axs[2,0].legend()

# set axis limits
axs[0,0].set_xlim(-1.5, 1.5)
axs[1,0].set_xlim(-1.5, 1.5)
axs[2,0].set_xlim(-1.5, 1.5)

plt.show()