import numpy as np
import matplotlib.pyplot as plt
import math

files = ['../Quartic/t1.dat',  '../Quartic/t15.dat', '../Quartic/t2.dat', '../Quartic/t25.dat']

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

# for the first subplot
for i in range(len(files)):

#    t = -1.0 - 0.5*i
    t = [-1, -1.5, -2, -2.5][i]

    axs[i].hist(data[i], bins=100, range=(-2.5, 2.5), density=True, color='gray',edgecolor='black',linewidth=0.5)
    axs[i].set_title("t = " + str(t))

    if i < 2:
        # plot theoretical distribution
        bt = np.sqrt(1/3 * (- 2*t + 2*np.sqrt(t**2 + 12)))
        ct = np.sqrt(1/2 * bt**2 + t)

        #from -bt, bt
        x = np.linspace(-bt, bt, 1000)
        y = 1/(2*np.pi) * np.sqrt(bt**2 - x**2) * (x**2 + ct**2)

        axs[i].plot(x, y, color='red', alpha=0.5, label='Densidade Teórica')

    else:
        # plot theoretical distribution
        at = np.sqrt(-2-t)
        bt = np.sqrt(2-t)

        x1 = np.linspace(-bt, -at, 500)
        x2 = np.linspace(at, bt, 500)
        y1 = np.abs(x1)/(2*np.pi) * np.sqrt((bt**2 - x1**2)*(x1**2 - at**2))
        y2 = np.abs(x2)/(2*np.pi) * np.sqrt((bt**2 - x2**2)*(x2**2 - at**2))

        axs[i].plot(x1, y1, color='red', alpha=0.8, label='Densidade Teórica')
        axs[i].plot(x2, y2, color='red', alpha=0.8)

    # label
    axs[i].legend(loc = 'upper right')

axs[0].set_ylabel("Frequencia")
axs[1].set_ylabel("Frequencia")
axs[2].set_ylabel("Frequencia")
axs[3].set_ylabel("Frequencia")
axs[len(files)-1].set_xlabel("Posição")

plt.show()