import numpy as np
import matplotlib.pyplot as plt

# Constants
cube_length_m = 0.5  # length of the cube, m
cylinder_radius_m = 0.1*cube_length_m  # radius of the cylinder, m
T_edge = 25.0  # temperature at the edges, °C
Q = 1  # total heat generation rate, W
# sigma = 0.0005  # standard deviation of the heat source, m
PDMS_thermal_conductivity_WpmK = 0.2
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_mps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)
print(PDMS_thermal_diffusivity_mps)
alpha = 1.172e-20  # thermal diffusivity of steel, m^2/s
Nx = 50  # number of grid points in the x direction
Ny = 50  # number of grid points in the y direction
Nz = 50  # number of grid points in the z direction
dx = cube_length_m / (Nx - 1)  # grid spacing in the x direction, m
dy = cube_length_m / (Ny - 1)  # grid spacing in the y direction, m
dz = cube_length_m / (Nz - 1)  # grid spacing in the z direction, m
dt = 0.01  # time step, s
x = np.linspace(0, cube_length_m, Nx)
y = np.linspace(0, cube_length_m, Ny)
z = np.linspace(0, cube_length_m, Nz)
X, Y, Z = np.meshgrid(x, y, z)

def compute_temperature_distribution_3D_cylinder_unbound(output_times):
    # Initial condition
    T = np.full((Nx, Ny, Nz), T_edge)

    # Heat source term
    q = Q / (np.sqrt(2 * np.pi) * cylinder_radius_m) * np.exp(-((X - cylinder_radius_m/2)**2 + (Y - cylinder_radius_m/2)**2))

    # Time stepping loop with specific output times
    output_indices = [int(t / dt) for t in output_times]
    output_temperatures = []

    for n in range(max(output_indices) + 1):
        # Compute temperature at next time step
        T_new = T + dt * (alpha * (np.roll(T, -1, axis=0) - 2*T + np.roll(T, 1, axis=0)) / dx**2
                          + alpha * (np.roll(T, -1, axis=1) - 2*T + np.roll(T, 1, axis=1)) / dy**2
                          + alpha * (np.roll(T, -1, axis=2) - 2*T + np.roll(T, 1, axis=2)) / dz**2
                          + q)
        # Enforce boundary conditions
        T_new[0, :, :] = T_edge
        T_new[-1, :, :] = T_edge
        T_new[:, 0, :] = T_edge
        T_new[:, -1, :] = T_edge
        T_new[:, :, 0] = T_edge
        # T_new[:, :, -1] = T_edge  # Do not enforce boundary condition on top face
        T = T_new
        # Save output if at specified time
        if n in output_indices:
            output_temperatures.append(T.copy())

    return output_temperatures

def plot_temperature_distribution_slices(X, Y, Z, output_temperatures, output_times):
    fig = plt.figure(figsize=(12, 8))
    min_temp = np.min(output_temperatures[0])
    max_temp = np.max(output_temperatures[-1])
    for i, T_out in enumerate(output_temperatures):
        ax = fig.add_subplot(2, 3, i+1, projection='3d')
        c = ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, Nz//2], levels=100, cmap='hot', vmin=min_temp, vmax=max_temp, alpha=0.75)
        ax.contourf(X[:, Ny//2, :], Y[:, Ny//2, :], T_out[:, Ny//2, :], zdir='y', offset=Y[0, Ny//2, 0], levels=100, cmap='hot', vmin=min_temp, vmax=max_temp, alpha=0.75)
        ax.contourf(X[Nx//2, :, :], Y[Nx//2, :, :], T_out[Nx//2, :, :], zdir='x', offset=X[Nx//2, 0, 0], levels=100, cmap='hot', vmin=min_temp, vmax=max_temp, alpha=0.75)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_zlabel('z (m)')
        ax.set_xlim([0, cube_length_m])
        ax.set_ylim([0, cube_length_m])
        ax.set_zlim([0, cube_length_m])
        fig.colorbar(c, ax=ax, label='Temperature (°C)')
    plt.tight_layout()
    plt.show()

# Compute and plot temperature distribution
output_times = [0, 1, 10]
output_temperatures = compute_temperature_distribution_3D_cylinder_unbound(output_times)
plot_temperature_distribution_slices(X, Y, Z, output_temperatures, output_times)