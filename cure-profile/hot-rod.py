import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Initialize array with n points all equal to 1
n = 100
array = np.ones(n)

# Apply arbitrary temperature profile function
tension = 1.5
array = (array * (np.exp(-np.linspace(0, 20, n)) ** tension) + (np.sin(np.linspace(0, np.pi * 3, n)) * 0.5) + 0.5) * 1000

# Set up the figure and axis
fig, ax = plt.subplots()
# Set the line color to blue and double the line width
line, = ax.plot(array, color='blue', linewidth=1.5) 
ax.set_title('Temperature evolution across a 1D rod')
ax.set_xlabel('l (cm)')
ax.set_ylabel('T (K)')

# Update function for animation
def update(i):
    global array

    # Roll laplacian
    array = (array + np.roll(array, 1) + np.roll(array, -1)) / 3

    # Adiabatic boundaries
    array[0] = (array[0] + array[1]) / 2
    array[-1] = (array[-1] + array[-2]) / 2

    # Update plot
    line.set_ydata(array)
    return line,

# Number of frames
t = 100

# Create animation
ani = FuncAnimation(fig, update, frames=t, blit=True)

# Save animation as a gif
ani.save('cure-profile\temperature_evolution.gif', writer='pillow', fps=20)

# Show the plot
plt.show()