import matplotlib.pyplot as plt
import numpy as np
from TracyWidom import TracyWidom
import math

files = ['../Tracy/sNm53.dat',
         '../Tracy/sNm43.dat',
         '../Tracy/sNm33.dat',
         '../Tracy/sNm23.dat',
         '../Tracy/sNm13.dat'
        ]   
        #  '../Tracy/sNm03.dat',
        #  '../Tracy/sN03.dat',
        #  '../Tracy/sN13.dat',
        #  '../Tracy/sN23.dat']

pots = ["-N^5/3", "-N^4/3", "-N^1","-N^2/3", "-N^1/3"] 
        #"-N^0", "N^0", "N^1/3", "N^2/3"]

data = []
beta = 2
N = 100

for file in files:
    dataaux = []
    with open(file) as f:
        for line in f:
            aux = line.split(" ")
            # filter empty string and \n
            aux = list(filter(lambda x: x != '\n' and x != '', aux))
            # make a list of floats
            aux = list(map(float, aux))
            #get max value
            dataaux.append(max(aux))
    data.append(dataaux)

# were doing tree subplots in one figure
fig, axs = plt.subplots(len(files), 1, sharey=True, tight_layout=True)

# title for the whole figure
fig.suptitle('Ensembles Tracy-Widom - Densidade')

# for the first subplot
for i in range(len(files)):

    adata = [(j * 2 / math.sqrt(2*beta)  - 2) * N**(2/3) for j in data[i]]

    # if data is out of 99 percentil, remove it
    adata = [i for i in adata if i < np.percentile(adata, 99)]

    axs[i].hist(adata, bins=100, range=(min(adata), max(adata)), density=True, color='gray', edgecolor='black', linewidth=0.5)
    axs[i].set_title("s = " + f"{pots[i]}")

    dataminrange = min(adata) - 0.1
    datamaxrange = max(adata) + 0.1
    print(dataminrange, datamaxrange, np.min(adata), np.max(adata))
    axs[i].set_xlim(-7, 1)

    
    #plot theoretical Tracy-Widom
    x = np.linspace(-7, 1, 1000)
    tw = TracyWidom(beta)
    y = [tw.pdf(i) for i in x]
    axs[i].plot(x, y, color='red', linewidth=1)


plt.show()
