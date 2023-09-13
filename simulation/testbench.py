import numpy as np
import matplotlib.pyplot as plt

# Constants
alpha = 0.2  # thermal diffusivity
dx = 0.01  # space step (m)
dy = dx
dt = 0.001  # time step (s)
nx = int(1/dx)  # number of points in x direction
ny = int(1/dy)  # number of points in y direction
nt = int(5/dt)  # number of timesteps

T = np.full((nx, ny), 25)  # initialize temperature array

# Generation of the constant 50 W in the leftmost 10% of the grid
# Assuming the heat generation is uniform over this region and per unit volume
q = 50 / (0.1 * 1 * 1)  # W/m^3

# Main loop
for t in range(nt+1):
    
    # Backup current temperature field
    Tn = T.copy()
    
    # Compute the temperature field at next time step
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            Txx = (Tn[i+1, j] - 2*Tn[i, j] + Tn[i-1, j]) / dx**2
            Tyy = (Tn[i, j+1] - 2*Tn[i, j] + Tn[i, j-1]) / dy**2
            T[i, j] = Tn[i, j] + alpha * dt * (Txx + Tyy)
            
            # Heat source in the leftmost 10%
            if i < nx/10:
                T[i, j] += q*dt
    
    # Boundary conditions
    T[:, 0] = 25  # left boundary
    T[:, -1] = 25  # right boundary
    T[0, :] = 25  # bottom boundary
    T[-1, :] = 25  # top boundary
    
    # Visualizing at times 0, 1 and 5 s
    if t == 0 or abs(t*dt - 1) < 1e-5 or abs(t*dt - 5) < 1e-5:
        plt.imshow(T, cmap='hot', extent=[0, 1, 0, 1])
        plt.colorbar()
        plt.title(f"Time = {t*dt:.2f} s")
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.show()

