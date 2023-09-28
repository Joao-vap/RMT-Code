#read file dataX.txt and plot histogram
import matplotlib.pyplot as plt
import numpy as np
from TracyWidom import TracyWidom

data = []
menores = []
with open('b4.txt') as f:
    for line in f:
        aux = line.split(" ")
        # filter empty string and \n
        aux = list(filter(lambda x: x != '\n' and x != '', aux))
        # make a list of floats
        aux = list(map(float, aux))
        for i in aux:
            data.append(i)
#       from aux append smallest value to menores
        menores.append(min(aux))


# were doing two subplots in one figure
fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

# for the first subplot
axs[0].hist(data, bins=50, range=(-1.5, 1.5), density=True, color='gray')
axs[0].set_title("Beta - Hermite | beta = 4")
axs[0].set_xlabel("Value")
axs[0].set_ylabel("Frequency")

# for the second subplot were doiong the tracywidow with menores data
# normalize menores to 1
menores = np.array(menores) + 1
plot, bins = np.histogram(menores, bins=50, range=(-3.1, 3.1), density=True)
# get bin centers
bin_centers = 0.5*(bins[1:] + bins[:-1])
# normalize histogram to 1
plot = plot/np.sum(plot)
# plot histogram in axes
axs[1].plot(bin_centers, plot, color='gray')
# name the axes
axs[1].set_title("Tracy-Widom | beta = 4")
axs[1].set_xlabel("Value")
axs[1].set_ylabel("Frequency")


# We are going to compute the density of eigenvalues of a Random Matrix
# in the models of GOE, GUE and GSE

import numpy as np
import matplotlib.pyplot as plt

# We define the number of eigenvalues
N = 10

# We define the number of realizations
M = 100000

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

# # We plot the theoretical density
# from TracyWidom import TracyWidom
# plt.plot(x, TracyWidom(1).pdf(x), label='TW1', color='red')
# plt.plot(x, TracyWidom(2).pdf(x), label='TW2', color='blue')
# plt.plot(x, TracyWidom(4).pdf(x), label='TW4', color='green')

axs[0].plot(x, density_GSE, label='GSE by RM', color='red')

# set x limits
axs[0].set_xlim(-1.5, 1.5)

# plot legend
axs[0].legend()

plt.show()
