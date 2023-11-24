import numpy as np
import matplotlib.pyplot as plt

##### Generate Array #####
# initialize array with n points all equal to 1
n = 100
array = np.ones(n)

# apply arbitrary temperature profile function
tension = 1.5
array = (array * (np.exp(-np.linspace(0, 20, n)) ** tension) + (np.sin(np.linspace(0, np.pi * 3, n)) * 0.5) + 0.5) * 1000


##### laplacian #####
# average each point with its neighbors with np.roll and show the graph at each step in the same window
t = 100
for i in range(t):
    # roll laplacian
    array = (array + np.roll(array, 1) + np.roll(array, -1)) / 3

    # adiabatic boundaries
    array[0] = (array[0] + array[1]) / 2
    array[-1] = (array[-1] + array[-2]) / 2

    # animate plot
    colormap = plt.cm.get_cmap('coolwarm')
    color = colormap(i / t)
    x_label = 'x'
    plt.plot(array, color=color, linewidth=.75)
    plt.title('Temperature evolution across a 1D rod')
    plt.xlabel('l (cm)')
    plt.ylabel('T (K)')
    plt.pause(0.05)


##### Visualize #####
plt.show()