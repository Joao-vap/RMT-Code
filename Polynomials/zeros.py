# we need do program elementary simmetric polynomials

# --------------------------------------
#imports
import matplotlib.pyplot as plt
import numpy as np
from mpmath import mp
#---------------------------------------

#first we need to import the data for the roots from file

filename = "./cilindrico/a5t1.dat"
fileroot = "./a5t1roots.dat"
N = 100
data = []

# with open(filename) as f:
#     for line in f:
#         aux = line.split(" ")
#         # filter empty string and \n
#         aux = list(filter(lambda x: x != '\n' and x != '', aux))
#         # make a list of complex numbers
#         aux = np.array([float(i) for i in aux])
#         #change shape to (N, 2)
#         aux = aux.reshape(-1, N)
#         dataaux = []
#         for i in range(aux.shape[1]):
#             dataaux.append(complex(aux[0, i], aux[1, i]))
#         data.append(dataaux)

# degree = len(data[0])

# def elementary_symmetric(vars, k):
#     n = len(vars)
#     dp = np.zeros((k+1, n+1), dtype=complex)
#     for j in range(n+1):
#         dp[0][j] = 1.0
#     for i in range(1, k+1):
#         for j in range(1, n+1):
#             dp[i][j] = dp[i][j-1] + vars[j-1] * dp[i-1][j-1]
#     return dp[k][n]

# coefs = [0 for _ in range(degree+1)]
# trials = 3200

# for j in range(trials):
#     listcoefs = np.poly(data[j])
#     for i in range(degree+1):
#         coefs[i] += listcoefs[i]
# #       coefs[i] += elementary_symmetric(data[j], i)
    
# # average
# for i in range(degree+1):
#     coefs[i] /= trials

# # plot 3d map of the polynomial
# x = np.linspace(-2, 2, 200)
# y = np.linspace(-2, 2, 200)
# X, Y = np.meshgrid(x, y)
# Z = np.zeros(X.shape)
# for i in range(X.shape[0]):
#     for j in range(X.shape[1]):
#         v = np.polyval(coefs, complex(X[i, j], Y[i, j])).real + np.polyval(coefs, complex(X[i, j], Y[i, j])).imag
#         Z[i, j] = v
    
# fig = plt.figure()
# # plot as colormap 2d countour
# plt.contour(X, Y, Z, [0, 0.5, 1, 10])
# #equal axis
# plt.axis('equal')
# plt.show()

# # find roots of polynomial x^degree - coefs[1] x^(degree-1) + coefs[2] x^(degree-2) - ... + (-1)^degree coefs[degree]
# # only real roots
# roots = np.roots(coefs)

# # plot roots
# plt.scatter([r.real for r in roots], [r.imag for r in roots])
# # equal axis
# plt.axis('equal')
# plt.show()

# # save roots in file
# with open(fileroot, 'w') as f:
#     for i in range(len(roots)):
#         f.write(str(roots[i].real) + " " + str(roots[i].imag) + "\n")


files = ["./a1t0roots.dat", "./a1t05roots.dat", "./a1t1roots.dat", "./a1t2roots.dat",
        "./a2t0roots.dat", "./a2t05roots.dat", "./a2t1roots.dat", "./a2t2roots.dat",
        "./a3t0roots.dat", "./a3t05roots.dat", "./a3t1roots.dat", "./a3t2roots.dat",
        "./a5t0roots.dat", "./a5t05roots.dat", "./a5t1roots.dat", "./a5t2roots.dat"]

measure = ["./cilindrico/a1t0.dat", "./cilindrico/a1t05.dat", "./cilindrico/a1t1.dat", "./cilindrico/a1t2.dat",
           "./cilindrico/a2t0.dat", "./cilindrico/a2t05.dat", "./cilindrico/a2t1.dat", "./cilindrico/a2t2.dat",
           "./cilindrico/a3t0.dat", "./cilindrico/a3t05.dat", "./cilindrico/a3t1.dat", "./cilindrico/a3t2.dat",
           "./cilindrico/a5t0.dat", "./cilindrico/a5t05.dat", "./cilindrico/a5t1.dat", "./cilindrico/a5t2.dat"]

#make plot 2x4
fig, axs = plt.subplots(4, 4, figsize=(15, 15))
fig.suptitle('Roots of the polynomials', fontsize=16)

for k in range(4):
    for j in range(4):

        rootname = files[k*4 + j]
        measurename = measure[k*4 + j]
        dataroots = []
        datameasure = []    

        with open(rootname) as r:
            for line in r:
                aux = line.split(" ")
                # filter empty string and \n
                aux = list(filter(lambda x: x != '\n' and x != '', aux))
                # make a list of complex numbers
                aux = np.array([float(i) for i in aux])
                # change shape to (N, 2)
                dataaux = complex(aux[0], aux[1])
                dataroots.append(dataaux)
        
        with open(measurename) as m:
            for line in m:
                aux = line.split(" ")
                # filter empty string and \n
                aux = list(filter(lambda x: x != '\n' and x != '', aux))
                # make a list of complex numbers
                aux = np.array([float(i) for i in aux])
                # change shape to (N, 2)
                aux = aux.reshape(-1, 100)
                for i in range(aux.shape[1]):
                    datameasure.append(complex(aux[0, i], aux[1, i]))

        trials = int(len(datameasure) / 100)
        degree = int(len(datameasure) / trials)
        coefs = [0 for _ in range(degree+1)]

        for q in range(trials):
            listcoefs = np.poly(datameasure[q*100:(q+1)*100])
            for r in range(degree):
                coefs[r] += listcoefs[r]
            
        # average
        for q in range(degree+1):
            coefs[q] /= trials
        
        minx = min([l.real for l in datameasure]) - 0.3
        maxx = max([l.real for l in datameasure]) + 0.3
        miny = min([l.imag for l in datameasure]) - 0.3
        maxy = max([l.imag for l in datameasure]) + 0.3

        #plot roots
        axs[k, j].scatter([l.real for l in dataroots], [l.imag for l in dataroots], color='red', s=2)
        #set axis limi and proportion
        axs[k, j].set_xlim(minx, maxx)
        axs[k, j].set_ylim(miny, maxy)
        axs[k, j].set_aspect('equal')

        s = 1500

        Z = np.zeros((s, s))
        for m in range(len(datameasure)):
            part = datameasure[m]
            coordx = int((part.real - minx) / (maxx - minx) * s)
            coordy = int((part.imag - miny) / (maxy - miny) * s)
            Z[coordy, coordx] += 1

        # cut = 1 all that is bigger than 1
        Z = np.clip(Z, 0, 1)

        #Z = Z / np.sum(Z)
        # with low alpha
        axs[k, j].imshow(1-Z, cmap='gray', extent=[minx, maxx, miny, maxy], alpha=0.3)

        # plot 3d map of the polynomial
        x = np.linspace(-2, 2, 200)
        y = np.linspace(-2, 2, 200)
        X, Y = np.meshgrid(x, y)
        Z = np.zeros(X.shape)
        for a in range(X.shape[0]):
            for b in range(X.shape[1]):
                pv = np.polyval(coefs, complex(X[a, b], Y[a, b]))
                v = pv.real
                Z[a, b] = v

        axs[k,j].contour(X, Y, Z, [0, 0.1, 0.5, 1, 5, 10, 1000]) 

# ylabel
axs[0, 0].set_ylabel('a1')
axs[1, 0].set_ylabel('a2')
axs[2, 0].set_ylabel('a3')
axs[3, 0].set_ylabel('a5')

# xlabel
axs[3, 0].set_xlabel('t0')
axs[3, 1].set_xlabel('t05')
axs[3, 2].set_xlabel('t1')
axs[3, 3].set_xlabel('t2')

# save plot
plt.savefig('roots.png')
plt.show()