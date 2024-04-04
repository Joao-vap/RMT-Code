#read file dataX.txt and plot histogram
import matplotlib.pyplot as plt
import numpy as np
from TracyWidom import TracyWidom
import math

#########################################################################

# files = ['../Gaussian/b1N10.dat', '../Gaussian/b1N100.dat',
#          '../Gaussian/b2N10.dat', '../Gaussian/b2N100.dat',
#          '../Gaussian/b4N10.dat', '../Gaussian/b4N100.dat',]

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

#     for j in range(2):

#         if j == 0:
#             N = 10
#             adata = [i / math.sqrt(2*beta) for i in data[2*i + j]]
#         else:
#             N = 100
#             adata = [i / math.sqrt(2*beta) for i in data[2*i + j]]

#         axs[i, j].hist(adata, bins=50, range=(-1.5, 1.5), density=True, color='gray', edgecolor='black',linewidth=0.5)
#         axs[i, j].set_title("beta = " + str(beta) + ", N = " + str(N))
#         axs[i, j].set_xlim(-1.5, 1.5)

# files = ['../Gaussian/b1N100.dat',
#          '../Gaussian/b2N100.dat',
#          '../Gaussian/b4N100.dat',]

# data = []

# N = 100

# for file in files:
#     dataaux = []
#     with open(file) as f:
#         for line in f:
#             aux = line.split(" ")
#             # filter empty string and \n
#             aux = list(filter(lambda x: x != '\n' and x != '', aux))
#             # make a list of floats
#             aux = list(map(float, aux))
#             # for i in aux:
#             #     dataaux.append(i)
#             # append max value in the list
#             dataaux.append(max(aux))
#     data.append(dataaux)

# for i in range(3):

#     if i == 0:
#         beta = 1
#     elif i == 1:
#         beta = 2
#     else:
#         beta = 4

# #    adata = [(j - math.sqrt(2*beta)) * N**(4/6) for j in data[i]]
#     adata = [(j * 2 / math.sqrt(2*beta)  - 2) * N**(2/3) for j in data[i]]

#     # center the data
#     #eliminate outliers
#     adata = [i for i in adata if i < 5 and i > -5]
#     adata = [i - np.mean(adata) for i in adata]

#     axs[i, 2].hist(adata, bins=80, range=(-5, 5), color='gray', density=True, edgecolor='black',linewidth=0.5)
#     axs[i, 2].set_title("Maior autovalor - beta = " + str(beta) + f", N = {N}")
#     axs[i, 2].set_xlim(-4, 4)


# axs[1, 0].set_ylabel("Frequência")
# axs[2, 1].set_xlabel("Posição")

# # plot the semi-circle for the axs[0,2], axs[1,2] and axs[2,2]
# x = np.linspace(-3, 3, 1000)
# y = 2/np.pi * np.sqrt(1 - x**2)

# axs[0, 1].plot(x, y, color='red', alpha=0.8,label='Wigner Semi-Circle')
# axs[1, 1].plot(x, y, color='red',alpha=0.8, label='Wigner Semi-Circle')
# axs[2, 1].plot(x, y, color='red', alpha=0.8,label='Wigner Semi-Circle')

# axs[0, 1].legend()
# axs[1, 1].legend()
# axs[2, 1].legend()

# x = np.linspace(-4, 4, 1000)

# # # We plot the theoretical density centered
# y1 = TracyWidom(1).pdf(x - 1.2065335745820)
# y2 = TracyWidom(2).pdf(x - 1.771086807411)
# y4 = TracyWidom(4).pdf(x - 2.306884893241)
# # y1 = TracyWidom(1).pdf(x)
# # y2 = TracyWidom(2).pdf(x)
# # y4 = TracyWidom(4).pdf(x)
# axs[0, 2].plot(x, y1, label='Tracy-Widow', color='green', alpha=0.8)
# axs[1, 2].plot(x, y2, label='Tracy-Widow', color='green', alpha=0.8)
# axs[2, 2].plot(x, y4, label='Tracy-Widow', color='green', alpha=0.8)


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

# files = ['../Quartic/t1.dat', '../Quartic/t2.dat', '../Quartic/t3.dat']

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

# # were doing two subplots in one figure
# fig, axs = plt.subplots(3, 2, sharey=True, tight_layout=True)

# # title for the whole figure
# fig.suptitle('Potenciais - Densidades')
# axs[0][0].set_title('Potencial Quártico')
# axs[1][0].set_title('Potencial Mônico')


# # for the first subplot
# for i in range(len(files)):

# #    t = -1.0 - 0.5*i
#     t = [-1, -2, -3][i]

#     axs[i][0].hist(data[i], bins=50, range=(-2.5, 2.5), density=True, color='gray',edgecolor='black',linewidth=0.5)
#     axs[i][0].set_title("t = " + str(t))

#     if i < 2:
#         # plot theoretical distribution
#         bt = np.sqrt(1/3 * (- 2*t + 2*np.sqrt(t**2 + 12)))
#         ct = np.sqrt(1/2 * bt**2 + t)

#         #from -bt, bt
#         x = np.linspace(-bt, bt, 1000)
#         y = 1/(2*np.pi) * np.sqrt(bt**2 - x**2) * (x**2 + ct**2)

#         axs[i][0].plot(x, y, color='red', alpha=0.5, label='Densidade Teórica')

#     else:
#         # plot theoretical distribution
#         at = np.sqrt(-2-t)
#         bt = np.sqrt(2-t)

#         x1 = np.linspace(-bt, -at, 500)
#         x2 = np.linspace(at, bt, 500)
#         y1 = np.abs(x1)/(2*np.pi) * np.sqrt((bt**2 - x1**2)*(x1**2 - at**2))
#         y2 = np.abs(x2)/(2*np.pi) * np.sqrt((bt**2 - x2**2)*(x2**2 - at**2))

#         axs[i][0].plot(x1, y1, color='red', alpha=0.8, label='Densidade Teórica')
#         axs[i][0].plot(x2, y2, color='red', alpha=0.8)

#     # label
#     axs[i][0].legend(loc = 'upper right')

# axs[1][0].set_ylabel("Frequencia")
# axs[len(files)-1][0].set_xlabel("Posição")


# # #########################################################################

# # #########################################################################

# files = ['../Monic/a1.dat', '../Monic/a3.dat', '../Monic/a5.dat']

# beta = 2
# t = 1

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

# def h(x, alpha, a):

#     y = x**(2*alpha - 2)

#     for j in range(1,alpha):

#         prod = 1
#         for l in range(1,j+1):
#             prod *= (2*l - 1)/(2*l)

#         y += x**(2*alpha-2-2*j) * a**(2*j) * prod

#     return y

# # for the first subplot
# for i in range(len(files)):

#     alpha = 2*i + 1
#     axs[i][1].hist(data[i], bins=50, range=(-2.5, 2.5), density=True, color='gray',edgecolor='black',linewidth=0.5)
#     axs[i][1].set_title("m = " + str(alpha))

#     # plot horizontal line in a and -a
#     a = alpha
#     for j in range(1,alpha+1):
#         a = a * (2*j-1)/(2*j)
#     a = a**(-1/(2*alpha))

#     # plot theoretical
#     x = np.linspace(-a, a, 1000)
#     y = [0 for i in range(len(x))]
#     for value in range(len(x)):
#         if x[value]**2 - a**2 < 0:
#             y[value] = alpha/math.pi * math.sqrt(a**2 - x[value]**2) * h(x[value], alpha, a)

#     axs[i][1].plot(x, y, color='red', alpha=0.8, label='Densidade Teórica')
#     axs[i][1].legend(loc = 'upper right')

# axs[len(files)-1][1].set_xlabel("Posição")


#########################################################################

# files = ['../V3/a2t0.dat', '../V3/a0t0.dat', '../V3/am2t0.dat',
#          '../V3/am2tm4.dat', '../V3/a0tm4.dat', '../V3/a2tm4.dat',]

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

# # were doing two subplots in one figure
# fig, axs = plt.subplots(2, 3, sharey=True, tight_layout=True)

# # title for the whole figure
# fig.suptitle('Densidades')

# # for the first subplot
# for i in range(2):

#     for j in range(3):

#         q = -2 + j*2
#         p = 0 - i*4

#         axs[i][j].hist(data[3*i + j], bins=50, range=(-2.5, 2.5), density=True, color='gray',edgecolor='black',linewidth=0.5)
#         axs[i][j].set_title(f"p = {p}, q = {q}")

# axs[0][0].set_ylabel("Frequencia")
# axs[1][0].set_ylabel("Frequencia")
# axs[1][0].set_xlabel("Posição")
# axs[1][1].set_xlabel("Posição")
# axs[1][2].set_xlabel("Posição")


#########################################################################

# # We are going to compute the density of eigenvalues of a Random Matrix
# # in the models of GOE, GUE and GSE

# # We define the number of eigenvalues
# N = 10

# # We define the number of realizations
# M = 4000

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

########################################################################

# # We plot the theoretical density
# from TracyWidom import TracyWidom
# plt.plot(x, TracyWidom(1).pdf(x), label='TW1', color='red')
# plt.plot(x, TracyWidom(2).pdf(x), label='TW2', color='blue')
# plt.plot(x, TracyWidom(4).pdf(x), label='TW4', color='green')

# ########################################################################

# # plor dataH.txt

# file = '../dataH.txt'

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

# # plot temporal series
# fig, ax = plt.subplots(1, 1, tight_layout=True)
# # scatter plot
# #ax.scatter(range(len(data)), data, color='blue', s=1)
# ax.plot(data, color='blue', alpha=0.2)
# # plor density
# #ax.hist(data, bins=50, range=(0.28, 0.32), density=True, color='gray')

# ax.set_title("Time Series")
# ax.set_xlabel("Time")
# ax.set_ylabel("Value")


# #########################################################################

# plot colormap of grid (NxN), N = 50
N = 100
data = []

filei = '../Complex/testi.dat'
filer = '../Complex/testr.dat'
file = '../Complex/test.dat'

datai = []
datar = []
dataV = []

with open(filei) as f:
    for line in f:
        aux = line.split(" ")
        # filter empty string and \n
        aux = list(filter(lambda x: x != '\n' and x != '', aux))
        # make a list of floats
        aux = list(map(float, aux))
        for i in aux:
            i = 1 if abs(i) > 1 else i
            datai.append(i)

with open(filer) as f:
    for line in f:
        aux = line.split(" ")
        # filter empty string and \n
        aux = list(filter(lambda x: x != '\n' and x != '', aux))
        # make a list of floats
        aux = list(map(float, aux))
        for i in aux:
            i = 1 if abs(i) > 1 else i
            datar.append(i)

with open(file) as f:
    for line in f:
        aux = line.split(" ")
        # filter empty string and \n
        aux = list(filter(lambda x: x != '\n' and x != '', aux))
        # make a list of floats
        aux = list(map(float, aux))
        for i in aux:
            i = 1 if i > 1 else i
            dataV.append(i)

# reshape data
dataaux = np.zeros((N,N))
for j in range(N):
    for i in range(N):
        iaux = N - i - 1
        jaux = j
        dataaux[iaux][jaux] = datai[N*j + i]
datai = dataaux

dataaux = np.zeros((N,N))
for j in range(N):
    for i in range(N):
        iaux = N - i - 1
        jaux = j
        dataaux[iaux][jaux] = datar[N*j + i]
datar = dataaux

dataaux = np.zeros((N,N))
for j in range(N):
    for i in range(N):
        iaux = N - i - 1
        jaux = j
        dataaux[iaux][jaux] = dataV[N*j + i]
dataV = dataaux

# plot colormap
fig, ax = plt.subplots(1, 3, tight_layout=True)

#set proportion
plt.gca().set_aspect('equal')

# plot real part
cax = ax[0].imshow(datar, cmap='hot', interpolation='nearest')
fig.colorbar(cax, ax=ax[0])

# plot imaginary part
cax = ax[1].imshow(datai, cmap='hot', interpolation='nearest')
fig.colorbar(cax, ax=ax[1])

# plot potential
cax = ax[2].imshow(dataV, cmap='hot', interpolation='nearest')
fig.colorbar(cax, ax=ax[2])


# data = np.zeros((10,10), dtype=complex)

# for i in range(10):
#     for j in range(10):
#         data[i][j] = complex(datar[i][j], datai[i][j])

# #plot vector field
# fig, ax = plt.subplots(1, 1, tight_layout=True)
# #set proportion
# plt.gca().set_aspect('equal')

# x = np.linspace(0, N, 1000)
# y = 0.2243545*(x-(N-1)/2) - (N-1)/2 
# ax.plot(x, y, color='red')


# for i in range(N):
#     for j in range(N):
#         ax.quiver(j, -i, datar[i][j], datai[i][j], color='blue', width=0.001)

# # #vector field
# # ax.quiver(datar, -datai, color='blue')

# # ax.set_title("Complex Variables")


########################################################################

#########################################################################

#plot two dimensional data

# files = ['../debug.dat']

# beta = 2

# data = []

# for file in files:
#     dataaux = []
#     with open(file) as f:
#         for line in f:
#             aux = line.split(" ")
#             # filter empty string and \n
#             aux = list(filter(lambda x: x != '\n' and x != '', aux))
#             # make a list of complex numbers
#             aux = [float(i) for i in aux]
#             for i in range(0, len(aux), 2):
#                 dataaux.append(complex(aux[i], aux[i+1]))
#     data.append(dataaux)

# # plot two dimensional data
# fig, ax = plt.subplots(1, 1, tight_layout=True)
# # set proportion
# plt.gca().set_aspect('equal')

# # scatter
# ax.scatter([i.real for i in data[0]], [i.imag for i in data[0]], color='grey', s=1)

# # plot unit circle

# # plot density lines
# # x = np.linspace(0, 2*np.pi, 1000)
# # y = np.linspace(0, 2*np.pi, 1000)

# # each = 1
# # for radius in range(0, 3000):
# #     print(radius)
# #     r = radius/1000
# #     count how many points are inside the circle and plot every 10%
# #     count = len([i for i in data[0] if abs(i) < r])
# #     if count > 0.25*each*len(data[0]):
# #         print(count, len(data[0]), count/len(data[0]))
# #         each += 1
# #         ax.plot(r*np.cos(x), r*np.sin(y), color='red')

# #     count = 0

# # ax.set_title("Complex Variables")


# #histogram of real part
# ax.hist([i.real for i in data[0]], bins=50, range=(-2.5, 2.5), density=True, color='red', edgecolor='black',linewidth=0.5, alpha=0.3)

# #histogram of real part
# ax.hist([i.imag for i in data[0]], bins=50, range=(-2.5, 2.5), density=True, color='blue', edgecolor='black',linewidth=0.5, alpha=0.3)

###########################################################################

# # extract data from energy
# # mean, median, std, max, min, skew, kurtosis

# files = ['../LogZnb/Hb1N50.dat', '../LogZnb/Hb1N100.dat', '../LogZnb/Hb1N150.dat',
#          '../LogZnb/Hb2N50.dat', '../LogZnb/Hb2N100.dat', '../LogZnb/Hb2N150.dat',
#          '../LogZnb/Hb4N50.dat', '../LogZnb/Hb4N100.dat', '../LogZnb/Hb4N150.dat']

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

# # plot energy data
# fig, ax = plt.subplots(1, 1, tight_layout=True)

# # scatter plot
# ax.scatter(range(len(data[7])), data[7], color='blue', s=1)

# ax.set_title("Energy")
# ax.set_xlabel("Time")
# ax.set_ylabel("Value")

# # extract mean, median, std, max, min, skew, kurtosis
# cutoff = [0.45, 0.45, 0.45]
# for i in range(len(files)):

#     # exclude outliers
#     data[i] = [j for j in data[i] if j < cutoff[i // 4] and j > -cutoff[i // 4]]

#     mean = np.mean(data[i])
#     median = np.median(data[i])
#     std = np.std(data[i])
#     max = np.max(data[i])
#     min = np.min(data[i])
#     skew = np.mean((data[i] - mean)**3) / std**3
#     kurtosis = np.mean((data[i] - mean)**4) / std**4

#     print(f'''File {files[i]}:
#     Mean: {mean}
#     Median: {median}
#     Std: {std}
#     Max: {max}
#     Min: {min}
#     Skew: {skew}
#     Kurtosis: {kurtosis}
#     -------------------------------------------------
#     ''')


###########################################################################

plt.show()
