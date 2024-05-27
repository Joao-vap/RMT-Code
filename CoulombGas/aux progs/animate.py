import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# we are going to animate many point dynamics on the line
# accordind to the file dataX.txt

datahist = []
data = []
with open('b2.txt') as f:
    l = []
    for line in f:
        aux = line.split(" ")
        # filter empty string and \n
        aux = list(filter(lambda x: x != '\n' and x != '', aux))
        for i in aux:
            l.append(float(i))
            datahist.append(float(i))
        data.append(l)
        l = []

# we need to create a figure and an axis
fig, ax = plt.subplots()

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown', 'pink']

# we need to create a function that will be called
# at each frame
def animate(i):
    # we need to clear the axis
    plt.cla()
    # we need to plot the points and color each with a different color
    for j in range(len(data[i])):
        ax.scatter(data[i][j], -0.1, color=colors[j % len(colors)])

    #plot histogram
    ax.hist(datahist[:i], bins=50, range=(-1.1, 1.1), density=True, color='gray')

    # we need to set the axis limits
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-0.4, 0.8)
    # we need to set the title
    ax.set_title("Beta - Hermite | beta = 2")
    # we need to set the labels
    ax.set_xlabel("Value")


# we need to create the animation
# we need to specify the figure, the function to call
# at each frame, the number of frames and the interval
# between frames in milliseconds
ani = animation.FuncAnimation(fig, animate, frames=len(data), interval=2000)

# we need to save the animation as mp4 video file
ani.save('b2.mp4', writer='ffmpeg', fps=10)

