#read file dataX.txt and plot histogram
import matplotlib.pyplot as plt
import numpy as np

data = []
with open('dataX.txt') as f:
    for line in f:
        aux = line.split(" ")
        # filter empty string and \n
        aux = list(filter(lambda x: x != '\n' and x != '', aux))
        for i in aux:
            data.append(float(i))

# histogram from -3 to 3 with 100 bins
plt.hist(data, bins=100, range=(-1.1, 1.1), density=True)
plt.title("Beta - Hermite | beta = 1")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()





