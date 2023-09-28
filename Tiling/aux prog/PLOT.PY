# PLOT THE FILE diamond.txt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Read the file
data = np.loadtxt('diamond.txt', skiprows=1)

m = '|'

t1 = mpl.markers.MarkerStyle(marker=m)
t1._transform = t1.get_transform().rotate_deg(45)

t2 = mpl.markers.MarkerStyle(marker=m)
t2._transform = t2.get_transform().rotate_deg(-45)

# Plot the data as an image with a colorbar with colors blue, green, yellow, red
for i in range(len(data)):
    for j in range(len(data[i])):
        if data[i][j] == 0:
            continue
        elif data[i][j] == 1:
            plt.plot(j, i, marker=t1, color='blue')
        elif data[i][j] == -1:
            plt.plot(j, i, marker=t1, color='green')
        elif data[i][j] == 2:
            plt.plot(j, i, marker=t2, color='red')
        elif data[i][j] == -2:
            plt.plot(j, i, marker=t2, color='black')

# plot circle with center in (len(data)/2, len(data)/2) and radius len(data)/2
circle = plt.Circle((len(data)/2, len(data)/2), len(data)/2, color='black', fill=False)
plt.gcf().gca().add_artist(circle)

# Add a title
plt.title('Diamond')

# set aspect ratio to 1
plt.axis('equal')

# set size o image
plt.gcf().set_size_inches(70, 70)

# Show the plot
plt.savefig('ad_circle.png', dpi=300)