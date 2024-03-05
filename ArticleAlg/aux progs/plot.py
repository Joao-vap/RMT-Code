#read file dataX.txt and plot histogram
import matplotlib.pyplot as plt
import numpy as np
from TracyWidom import TracyWidom

#########################################################################

# files = ['../Ensemble Data/b1N10.txt', '../Ensemble Data/b1N50.txt', '../Ensemble Data/b1N100.txt',
#          '../Ensemble Data/b2N10.txt', '../Ensemble Data/b2N50.txt', '../Ensemble Data/b2N100.txt',
#          '../GSE/n10.dat', '../Ensemble Data/b4N50.txt', '../Ensemble Data/b4N100.txt']

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
# fig.suptitle('Gaussian Ensembles - Density')

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
#             N = 50
#         else:
#             N = 100

#         axs[i, j].hist(data[3*i + j], bins=50, range=(-1.5, 1.5), density=True, color='gray')
#         axs[i, j].set_title("beta = " + str(beta) + ", N = " + str(N))
#         axs[i, j].set_xlabel("Value")
#         axs[i, j].set_ylabel("Frequency")
#         axs[i, j].set_xlim(-1.5, 1.5)
        
# # plot the semi-circle for the axs[0,2], axs[1,2] and axs[2,2]
# x = np.linspace(-1.5, 1.5, 1000)
# y = 2/np.pi * np.sqrt(1 - x**2)

# axs[0, 2].plot(x, y, color='red')
# axs[1, 2].plot(x, y, color='red')
# axs[2, 2].plot(x, y, color='red')

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
fig, axs = plt.subplots(len(files), 1, sharey=True, tight_layout=True)

# title for the whole figure
fig.suptitle('Quartic and Monic Potential - Density (N=50)')

# for the first subplot
for i in range(len(files)):

#    t = -1.0 - 0.5*i
    t = [-1, -1.5, -2, -2.5, -3][i]

    axs[i].hist(data[i], bins=100, range=(-2.5, 2.5), density=True, color='gray')
    axs[i].set_title("t = " + str(t))
    axs[i].set_ylabel("Frequency")

    if i < 2:
        # plot theoretical distribution
        bt = np.sqrt(1/3 * (- 2*t + 2*np.sqrt(t**2 + 12)))
        ct = np.sqrt(1/2 * bt**2 + t)

        #from -bt, bt
        x = np.linspace(-bt, bt, 1000)
        y = 1/(2*np.pi) * np.sqrt(bt**2 - x**2) * (x**2 + ct**2)

        axs[i].plot(x, y, color='red', alpha=0.5)

    else:
        # plot theoretical distribution
        at = np.sqrt(-2-t)
        bt = np.sqrt(2-t)

        x1 = np.linspace(-bt, -at, 500)
        x2 = np.linspace(at, bt, 500)
        y1 = np.abs(x1)/(2*np.pi) * np.sqrt((bt**2 - x1**2)*(x1**2 - at**2))
        y2 = np.abs(x2)/(2*np.pi) * np.sqrt((bt**2 - x2**2)*(x2**2 - at**2))

        axs[i].plot(x1, y1, color='red', alpha=0.5)
        axs[i].plot(x2, y2, color='red', alpha=0.5)


axs[len(files)-1].set_xlabel("Value")


#########################################################################

#########################################################################

# files = ['../Monic/M1.txt', '../Monic/M2.txt', '../Monic/M3.txt',
#             '../Monic/M4.txt', '../Monic/M5.txt']

# beta = 2

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

# # # were doing two subplots in one figure
# # fig, axs = plt.subplots(5, 1, sharey=True, tight_layout=True)

# # # title for the whole figure
# # fig.suptitle('Monic Potential - Density (N=100)')

# # for the first subplot
# for i in range(5):

#     alpha = i + 1
#     axs[i][1].hist(data[i], bins=500, range=(-2.5, 2.5), density=True, color='gray')
#     axs[i][1].set_title("alpha = " + str(alpha))

#     # plot horizontal line in a and -a
#     a = 0.5
#     for j in range(1,alpha):
#         a = a * (2*j-1)/(2*j)
#     a = a**(-1/(alpha))

#     axs[i][1].axvline(x=a, color='red', alpha=0.5)
#     axs[i][1].axvline(x=-a, color='red', alpha=0.5)

# axs[4][1].set_xlabel("Value")

# # plot the semi-circle in axs[0]
# x = np.linspace(-2, 2, 1000)
# y = 1/(2*np.pi) * np.sqrt(4 - x**2)

# axs[0][1].plot(x, y, color='blue', alpha=0.5)

#########################################################################

#########################################################################

# # We are going to compute the density of eigenvalues of a Random Matrix
# # in the models of GOE, GUE and GSE

# # We define the number of eigenvalues
# N = 10

# # We define the number of realizations
# M = 100000

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


# axs[0,0].plot(x, density_GOE, color='blue', alpha=0.4)
# axs[1,0].plot(x, density_GUE, color='blue', alpha=0.4)
# axs[2,0].plot(x, density_GSE, color='blue', alpha=0.4)

# # set axis limits
# axs[0,0].set_xlim(-3, 3)
# axs[1,0].set_xlim(-3, 3)
# axs[2,0].set_xlim(-3, 3)

#########################################################################

plt.show()
