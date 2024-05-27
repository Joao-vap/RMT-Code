# plot complex roots
import matplotlib.pyplot as plt
import numpy as np

# files = ["./a1t2roots.dat", "./a2t2roots.dat", "./a3t2roots.dat", "./a5t2roots.dat"]
# measure = ["./a1t2.dat", "./a2t2.dat", "./a3t2.dat", "./a5t2.dat"]

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

    