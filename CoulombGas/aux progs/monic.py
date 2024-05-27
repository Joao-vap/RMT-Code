import numpy as np
import matplotlib.pyplot as plt
import math

files = ['../Monic/a1.dat','../Monic/a2.dat', '../Monic/a3.dat', '../Monic/a4.dat']

fig, axs = plt.subplots(len(files), 1, sharey=True, tight_layout=True, figsize=(10, 10))

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
for i in range(len(files)):

    alpha = i + 1
    axs[i].hist(data[i], bins=100, range=(-2.5, 2.5), density=True, color='gray',edgecolor='black',linewidth=0.5)
    axs[i].set_title("m = " + str(alpha))

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

    axs[i].plot(x, y, color='red', alpha=0.8, label='Densidade Teórica')
    axs[i].legend(loc = 'upper right')

axs[len(files)-1].set_xlabel("Posição")
axs[0].set_ylabel("Frequencia")
axs[1].set_ylabel("Frequencia")
axs[2].set_ylabel("Frequencia")
axs[3].set_ylabel("Frequencia")

plt.show()
