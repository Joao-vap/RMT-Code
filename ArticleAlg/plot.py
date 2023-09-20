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
plt.hist(data, bins=100, range=(-2, 2), density=True)
plt.title("Histogram of dataX.txt")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()



