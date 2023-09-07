import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

## physical parameters ##
cube_length_m = 10  # Length of the cube, m
cylinder_radius_m = 0.05  # Radius of the cylinder, m
T_edge = 25.0  # Temperature at the edges, °C
T_air = 25.0  # Temperature of the surrounding air, °C
h_conv = 10.0  # Convective heat transfer coefficient, W/(m^2 K)
Q = 1  # Total heat generation rate, W
PDMS_thermal_conductivity_WpmK = 0.2e-10
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)
TD_steel = 1.172e-1  # Thermal diffusivity of steel, m^2/s

## mesh-grid parameters ##
Nx, Ny, Nz = 50, 50, 50  # Number of grid points in each dimension
dx = cube_length_m / (Nx - 1)  # Grid spacing in the x direction, m
dy = cube_length_m / (Ny - 1)  # Grid spacing in the y direction, m
dz = cube_length_m / (Nz - 1)  # Grid spacing in the z direction, m
dt = 0.05  # Time step, s
x = np.linspace(0, cube_length_m, Nx)
y = np.linspace(0, cube_length_m, Ny)
z = np.linspace(0, cube_length_m, Nz)
X, Y, Z = np.meshgrid(x, y, z)

## heat source with exponentially decaying power to emulate laser absorption ##
z_end_fraction = 0.2 # fraction of cube length where heat source is at 10% of original value
z_end = z_end_fraction * cube_length_m
abs_coeff = np.log(0.1) / -z_end
q = Q / (np.sqrt(2 * np.pi) * cylinder_radius_m) * np.exp(-((X - cube_length_m / 2)**2 + (Y - cube_length_m / 2)**2)) * np.exp(abs_coeff * (cube_length_m - Z))
q[Z < cube_length_m - z_end] = 0

def Compute_T_3D_RK(output_times):
    '''computes T at each grid point at each time step using Runge-Kutta method for input of times in seconds'''
    # Initialize temperature distribution
    T = np.full((Nx, Ny, Nz), T_edge)

    # Define heat source
    q = Q / (np.sqrt(2 * np.pi) * cylinder_radius_m) * np.exp(-((X - cube_length_m / 2)**2 + (Y - cube_length_m / 2)**2))

    # Time steps
    output_indices = [int(t / dt) for t in output_times]
    output_temperatures = []
    for n in range(max(output_indices) + 1):
        k1 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T, -1, axis=0) - 2*T + np.roll(T, 1, axis=0)) / dx**2
                   + PDMS_thermal_diffusivity_m2ps * (np.roll(T, -1, axis=1) - 2*T + np.roll(T, 1, axis=1)) / dy**2
                   + PDMS_thermal_diffusivity_m2ps * (np.roll(T, -1, axis=2) - 2*T + np.roll(T, 1, axis=2)) / dz**2
                   + q)
        
        T_temp = T + k1 / 2
        k2 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                   + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                   + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=2) - 2*T_temp + np.roll(T_temp, 1, axis=2)) / dz**2
                   + q)
        
        T_temp = T + k2 / 2
        k3 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                   + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                   + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=2) - 2*T_temp + np.roll(T_temp, 1, axis=2)) / dz**2
                   + q)
        
        T_temp = T + k3
        k4 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                   + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                   + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=2) - 2*T_temp + np.roll(T_temp, 1, axis=2)) / dz**2
                   + q)
        
        T_new = T + (k1 + 2*k2 + 2*k3 + k4) / 6
        print(T_new)
        # Apply boundary conditions with convection cooling at the top surface
        T_new[0, :, :] = T_edge
        T_new[-1, :, :] = T_edge
        T_new[:, 0, :] = T_edge
        T_new[:, -1, :] = T_edge
        T_new[:, :, 0] = T_edge
        T_new[:, :, -1] = T_new[:, :, -1] - h_conv * (T_new[:, :, -1] - T_air) * dt / (dz)
        T = T_new
        if n in output_indices:
            output_temperatures.append(T.copy())
    return output_temperatures

# Function to plot temperature slices
def Plot_T_Slices_RK(X, Y, Z, output_temperatures, output_times):
    fig = plt.figure(figsize=(12, 8))
    custom_cmap = plt.cm.get_cmap('hot', 20) # also, 'coolwarm'
    for i, T_out in enumerate(output_temperatures):
        ax = fig.add_subplot(1, len(output_times), i+1, projection='3d')
        c = ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, Nz//2], 
                        levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        # ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, Nz//2], 
        #                 levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        # ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, 0], 
        #                 levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        # ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='x', offset=0, 
        #                 levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_zlabel('z (m)')
        ax.set_xlim([0, cube_length_m])
        ax.set_ylim([0, cube_length_m])
        ax.set_zlim([0, cube_length_m])
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=np.arange(20, 116, 10))
    plt.tight_layout()
    plt.show()

def save_output_temperatures(output_temperatures, filename="output_temperatures.npz"):
    """
    Save the output temperatures to a file.
    
    Parameters:
    - output_temperatures: List of 3D numpy arrays containing temperature data.
    - filename: Name of the file to save the data to.
    """
    np.savez(filename, *output_temperatures)

output_times = [0, 5, 10, 30]
output_temperatures_RK = Compute_T_3D_RK(output_times)

## plot or save to file as npz ##
Plot_T_Slices_RK(X, Y, Z, output_temperatures_RK, output_times)
save_output_temperatures(output_temperatures_RK)