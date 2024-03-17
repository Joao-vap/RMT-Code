#read file dataX.txt and plot histogram
import matplotlib.pyplot as plt
import numpy as np
from TracyWidom import TracyWidom
import math

#########################################################################

# files = ['../Gaussian/b1N10.dat', '../Gaussian/b1N20.dat', '../Gaussian/b1N50.dat',
#          '../Gaussian/b2N10.dat', '../Gaussian/b2N20.dat', '../Gaussian/b2N50.dat',
#          '../Gaussian/b4N10.dat', '../Gaussian/b4N20.dat', '../Gaussian/b4N50.dat']

# data = []

# for file in files:
#     dataaux = []
#     with open(file) as f:
#         for line in f:
#             aux = line.split(" ")
#             # filter empty string and \n
#             aux = list(filter(lambda x: x != '\n' and x != '', aux))
#             # make a list of floats
#             aux = list(map(float, aux))
#             for i in aux:
#                 dataaux.append(i)
#     data.append(dataaux)

# # were doing two subplots in one figure
# fig, axs = plt.subplots(3, 3, sharey=True, tight_layout=True)

# # title for the whole figure
# fig.suptitle('Ensembles Gaussianos - Densidade')

# # for the first subplot
# for i in range(3):

#     if i == 0:
#         beta = 1
#     elif i == 1:
#         beta = 2
#     else:
#         beta = 4

#     for j in range(3):

#         if j == 0:
#             N = 10
#         elif j == 1:
#             N = 20
#         else:
#             N = 50

#         axs[i, j].hist(data[3*i + j], bins=50, range=(-1.5, 1.5), density=True, color='gray')
#         axs[i, j].set_title("beta = " + str(beta) + ", N = " + str(N))
#         axs[i, j].set_xlim(-1.5, 1.5)

# axs[1, 0].set_ylabel("Frequência")
# axs[2, 1].set_xlabel("Posição")
        
# # plot the semi-circle for the axs[0,2], axs[1,2] and axs[2,2]
# x = np.linspace(-1.5, 1.5, 1000)
# y = 2/np.pi * np.sqrt(1 - x**2)

# axs[0, 2].plot(x, y, color='red', alpha=0.8,label='Wigner Semi-Circle')
# axs[1, 2].plot(x, y, color='red',alpha=0.8, label='Wigner Semi-Circle')
# axs[2, 2].plot(x, y, color='red', alpha=0.8,label='Wigner Semi-Circle')

# # show labels
# axs[0, 2].legend()
# axs[1, 2].legend()
# axs[2, 2].legend()

#########################################################################

# # test gaussian variables

# file = '../gaussian.dat'

# data = []
# with open(file) as f:
#     for line in f:
#         aux = line.split(" ")
#         # filter empty string and \n
#         aux = list(filter(lambda x: x != '\n' and x != '', aux))
#         # make a list of floats
#         aux = list(map(float, aux))
#         for i in aux:
#             data.append(i)

# # plot density
# fig, ax = plt.subplots(1, 1, tight_layout=True)

# ax.hist(data, bins=50, range=(-3, 3), density=True, color='gray')
# ax.set_title("Gaussian Variables")
# ax.set_xlabel("Value")
# ax.set_ylabel("Frequency")
# ax.set_xlim(-3, 3)

# # plot teh gaussian density
# x = np.linspace(-3, 3, 1000)
# y = 1/np.sqrt(2*np.pi) * np.exp(-x**2/2)

# ax.plot(x, y, color='red')


#########################################################################

files = ['../Quartic/t1.dat', '../Quartic/t15.dat', '../Quartic/t2.dat', '../Quartic/t25.dat', '../Quartic/t3.dat']

beta = 2

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
fig, axs = plt.subplots(5, 2, sharey=True, tight_layout=True)

# title for the whole figure
fig.suptitle('Potencial Quártico e Mônico - Densidade (N=50)')

# for the first subplot
for i in range(len(files)):

#    t = -1.0 - 0.5*i
    t = [-1, -1.5, -2, -2.5, -3][i]

    axs[i][0].hist(data[i], bins=50, range=(-2.5, 2.5), density=True, color='gray')
    axs[i][0].set_title("t = " + str(t))

    if i < 2:
        # plot theoretical distribution
        bt = np.sqrt(1/3 * (- 2*t + 2*np.sqrt(t**2 + 12)))
        ct = np.sqrt(1/2 * bt**2 + t)

        #from -bt, bt
        x = np.linspace(-bt, bt, 1000)
        y = 1/(2*np.pi) * np.sqrt(bt**2 - x**2) * (x**2 + ct**2)

        axs[i][0].plot(x, y, color='red', alpha=0.5, label='Densidade Teórica')

    else:
        # plot theoretical distribution
        at = np.sqrt(-2-t)
        bt = np.sqrt(2-t)

        x1 = np.linspace(-bt, -at, 500)
        x2 = np.linspace(at, bt, 500)
        y1 = np.abs(x1)/(2*np.pi) * np.sqrt((bt**2 - x1**2)*(x1**2 - at**2))
        y2 = np.abs(x2)/(2*np.pi) * np.sqrt((bt**2 - x2**2)*(x2**2 - at**2))

        axs[i][0].plot(x1, y1, color='red', alpha=0.8, label='Densidade Teórica')
        axs[i][0].plot(x2, y2, color='red', alpha=0.8)

    # label
    axs[i][0].legend(loc = 'upper right')

axs[2][0].set_ylabel("Frequencia")
axs[len(files)-1][0].set_xlabel("Posição")


#########################################################################

#########################################################################

files = ['../Monic/a1.dat', '../Monic/a2.dat', '../Monic/a3.dat',
            '../Monic/a4.dat', '../Monic/a5.dat']

beta = 2
t = 1

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

# # were doing two subplots in one figure
# fig, axs = plt.subplots(5, 1, sharey=True, tight_layout=True)

# # title for the whole figure
# fig.suptitle('Monic Potential - Density (N=100)')

def h(x, alpha, a):

    y = x**(2*alpha - 2)

    for j in range(1,alpha):

        prod = 1
        for l in range(1,j+1):
            prod *= (2*l - 1)/(2*l)

        y += x**(2*alpha-2-2*j) * a**(2*j) * prod

    return y

# for the first subplot
for i in range(5):

    alpha = i + 1
    axs[i][1].hist(data[i], bins=50, range=(-2.5, 2.5), density=True, color='gray')
    axs[i][1].set_title("m = " + str(alpha))

    # plot horizontal line in a and -a
    a = alpha
    for j in range(1,alpha+1):
        a = a * (2*j-1)/(2*j)
    a = a**(-1/(2*alpha))

    # plot theoretical
    x = np.linspace(-a, a, 1000)
    y = [0 for i in range(len(x))]
    for value in range(len(x)):
        if x[value]**2 - a**2 < 0:
            y[value] = alpha/math.pi * math.sqrt(a**2 - x[value]**2) * h(x[value], alpha, a)

    axs[i][1].plot(x, y, color='red', alpha=0.8, label='Densidade Teórica')
    axs[i][1].legend(loc = 'upper right')

axs[4][1].set_xlabel("Posição")


#########################################################################

#########################################################################

# # We are going to compute the density of eigenvalues of a Random Matrix
# # in the models of GOE, GUE and GSE

# # We define the number of eigenvalues
# N = 10

# # We define the number of realizations
# M = 10000

# # Gaussian Orthogonal Ensemble
# def GOE(N, escale = False):
#     beta = 1
#     # We generate a random matrix
#     A = np.random.randn(N,N)
#     # We make it symmetric
#     A = (A + A.T)/2
#     # We compute the eigenvalues
#     eig = np.linalg.eigvalsh(A)
#     # return unique eigenvalues
#     if escale:
#         return np.unique(eig)/(np.sqrt(N*beta))
#     return np.unique(eig)

# # Gaussian Unitary Ensemble
# def GUE(N, escale = False):
#     beta = 2
#     # We generate a random matrix
#     A = np.random.randn(N,N) + 1j*np.random.randn(N,N)
#     # We make it Hermitian
#     A = (A + A.T.conj())/2
#     # We compute the eigenvalues
#     eig = np.linalg.eigvalsh(A)
#     # return unique eigenvalues
#     if escale:
#         return np.unique(eig)/(np.sqrt(N*beta))
#     return np.unique(eig)

# # Gaussian Symplectic Ensemble
# def GSE(N, escale = False):
#     beta = 4
#     # We generate a random matrix with quaternions entries
#     A = np.random.randn(N,N) + 1j*np.random.randn(N,N)
#     B = np.random.randn(N,N) + 1j*np.random.randn(N,N)
#     M = np.block([[A, -B.conj()], [B, A.conj()]])
#     # We make it Hermitian
#     M = (M + M.T.conj())/2
#     # We compute the eigenvalues
#     eig = np.linalg.eigvalsh(M)
#     # return unique eigenvalues
#     if escale:
#         return np.unique(eig)/(np.sqrt(N*beta))
#     return np.unique(eig)

# # We compute the eigenvalues of the random matrices
# eig_GOE = np.zeros((M,N))
# eig_GUE = np.zeros((M,N))
# eig_GSE = np.ones((M,2*N))*1000

# for i in range(M):
#     eig_GOE[i,:] = GOE(N, escale=True)
#     eig_GUE[i,:] = GUE(N, escale=True)
#     # complete GSE with 1000 as eigenvalues
#     value = GSE(N, escale=True)
#     eig_GSE[i,:len(value)] = value

# # We compute the density of eigenvalues
# x = np.linspace(-5,5,1000)
# dx = x[1] - x[0]

# density_GOE = np.zeros(x.shape)
# density_GUE = np.zeros(x.shape)
# density_GSE = np.zeros(x.shape)

# eig_GOE = eig_GOE / np.sqrt(2)
# eig_GUE = eig_GUE / np.sqrt(2)
# eig_GSE = eig_GSE / np.sqrt(2)

# for i in range(M):
#     density_GOE[:-1] += np.histogram(eig_GOE[i,:], bins=x)[0]
#     density_GUE[:-1] += np.histogram(eig_GUE[i,:], bins=x)[0]
#     density_GSE[:-1] += np.histogram(eig_GSE[i,:], bins=x)[0]

# density_GOE /= (M*N)/100
# density_GUE /= (M*N)/100
# density_GSE /= (M*N*2)/100

# # # We plot the density of eigenvalues
# # plt.plot(x, density_GOE, label='GOE', color='red', alpha=0.4)
# # plt.plot(x, density_GUE, label='GUE', color='blue', alpha=0.4)
# # plt.plot(x, density_GSE, label='GSE', color='green', alpha=0.4)

# # # We plot the theoretical density
# # from TracyWidom import TracyWidom
# # plt.plot(x, TracyWidom(1).pdf(x), label='TW1', color='red')
# # plt.plot(x, TracyWidom(2).pdf(x), label='TW2', color='blue')
# # plt.plot(x, TracyWidom(4).pdf(x), label='TW4', color='green')


# axs[0,0].plot(x, density_GOE, color='blue', alpha=0.4, label='Matrix Model')
# axs[1,0].plot(x, density_GUE, color='blue', alpha=0.4, label='Matrix Model')
# axs[2,0].plot(x, density_GSE, color='blue', alpha=0.4, label='Matrix Model')

# #legend
# axs[0,0].legend()
# axs[1,0].legend()
# axs[2,0].legend()

# # set axis limits
# axs[0,0].set_xlim(-1.5, 1.5)
# axs[1,0].set_xlim(-1.5, 1.5)
# axs[2,0].set_xlim(-1.5, 1.5)

#########################################################################

plt.show()
