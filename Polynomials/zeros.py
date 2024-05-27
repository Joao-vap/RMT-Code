# we need do program elementary simmetric polynomials

# --------------------------------------
#imports
import matplotlib.pyplot as plt
import numpy as np
from mpmath import mp
#---------------------------------------

#first we need to import the data for the roots from file

filename = "./cilindrico/a1t2.dat"
fileroot = "./a1t2roots.dat"
N = 100
data = []

with open(filename) as f:
    for line in f:
        aux = line.split(" ")
        # filter empty string and \n
        aux = list(filter(lambda x: x != '\n' and x != '', aux))
        # make a list of complex numbers
        aux = np.array([float(i) for i in aux])
        #change shape to (N, 2)
        aux = aux.reshape(-1, N)
        dataaux = []
        for i in range(aux.shape[1]):
            dataaux.append(complex(aux[0, i], aux[1, i]))
        data.append(dataaux)

degree = len(data[0])

def elementary_symmetric(vars, k):
    n = len(vars)
    dp = np.zeros((k+1, n+1), dtype=complex)
    for j in range(n+1):
        dp[0][j] = 1.0
    for i in range(1, k+1):
        for j in range(1, n+1):
            dp[i][j] = dp[i][j-1] + vars[j-1] * dp[i-1][j-1]
    return dp[k][n]

coefs = [0 for _ in range(degree+1)]
trials = 3200

for j in range(trials):
    listcoefs = np.poly(data[j])
    for i in range(degree+1):
        coefs[i] += listcoefs[i]
#       coefs[i] += elementary_symmetric(data[j], i)
    
# average
for i in range(degree+1):
    coefs[i] /= trials

# # plot 3d map of the polynomial
# x = np.linspace(-0.2, 0.2, 500)
# y = np.linspace(-0.2, 0.2, 500)
# X, Y = np.meshgrid(x, y)
# Z = np.zeros(X.shape)
# for i in range(X.shape[0]):
#     for j in range(X.shape[1]):
#         v = abs(np.polyval(coefs, complex(X[i, j], Y[i, j])))
#         Z[i, j] = v
    
# fig = plt.figure()
# # plot as colormap 2d countour
# plt.contour(X, Y, Z, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 1])
# #equal axis
# plt.axis('equal')
# plt.show()

# find roots of polynomial x^degree - coefs[1] x^(degree-1) + coefs[2] x^(degree-2) - ... + (-1)^degree coefs[degree]
# only real roots
roots = np.roots(coefs)

# plot roots
plt.scatter([r.real for r in roots], [r.imag for r in roots])
# equal axis
plt.axis('equal')
plt.show()

# save roots in file
with open(fileroot, 'w') as f:
    for i in range(len(roots)):
        f.write(str(roots[i].real) + " " + str(roots[i].imag) + "\n")