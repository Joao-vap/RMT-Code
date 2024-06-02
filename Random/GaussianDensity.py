# We are going to compute the density of eigenvalues of a Random Matrix
# in the models of GOE, GUE and GSE

import numpy as np
import matplotlib.pyplot as plt

# We define the number of eigenvalues
N = 8

# We define the number of realizations
M = 1000000

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

for i in range(M):
    density_GOE[:-1] += np.histogram(eig_GOE[i,:], bins=x)[0]
    density_GUE[:-1] += np.histogram(eig_GUE[i,:], bins=x)[0]
    density_GSE[:-1] += np.histogram(eig_GSE[i,:], bins=x)[0]

density_GOE /= (M*N)/100
density_GUE /= (M*N)/100
density_GSE /= (M*N*2)/100

# We plot the density of eigenvalues
plt.plot(x, density_GOE, label='GOE', color='red', alpha=0.4)
plt.plot(x, density_GUE, label='GUE', color='blue', alpha=0.4)
plt.plot(x, density_GSE, label='GSE', color='green', alpha=0.4)

# show x from -2 to 2
plt.xlim(-2,2)
# set grid
plt.grid()
plt.legend()
plt.show()
