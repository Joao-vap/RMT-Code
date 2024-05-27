# import data

import matplotlib.pyplot as plt
import numpy as np


with open('../Tracy/sN23.dat') as f:
    # read first line
    line = f.readline()
    print(line)
    # split line by space
    aux = line.split(" ")
    # filter empty string and \n
    aux = list(filter(lambda x: x != '\n' and x != '', aux))
    # make a list of floats
    aux = list(map(float, aux))

# plot scatter
plt.scatter(aux,[0 for i in range(len(aux))],color='gray', s=1)

plt.show()
    

