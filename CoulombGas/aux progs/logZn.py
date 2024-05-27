# import hamiltonian data from ../LogZnb/Hb2N100.dat

import numpy as np
import matplotlib.pyplot as plt
import math

# read data from file
file = '../LogZnb/Hb1N100.dat'
N = 100
beta = 1
with open(file, 'r') as f:
    data = f.readlines()
    #remove breakline
    data = [x.strip() for x in data]

#remove outliers
data = np.array([float(x) for x in data if float(x) < 1])

l = len(data)

# with first 2000 data
H1 = data[:math.floor(l/2)]
# with last 2000 data
H2 = data[-math.floor(l/2):]

# calculate logZn mean with logZn_med = N^2/S * sum(-beta * H1)
logZn = N**2 * np.mean(-beta*H1)
print('logZn mean:', logZn)

# calculate alpha mean with alpha_med = 1/S * sum((logZn + beta * H2 * N^2)/Nlog(N))
alpha = np.mean((logZn + beta*H2*N**2)/(N*math.log(N)))
print('alpha mean:', alpha)