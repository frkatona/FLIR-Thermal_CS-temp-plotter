import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from tqdm import tqdm

## Material physical parameters ##
square_length_m = 1  # m
T_edge = 25.0  # °C
T_air = 25.0  # °C
h_conv = 10.0  # W/(m^2 K)
Q = 50 # laser power, W  
PDMS_thermal_conductivity_WpmK = 0.2  # W/(m K)
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

## Mesh-grid parameters ##
Nx = 101  # Number of grid points in each dimension
Ny = Nx
dx = square_length_m / (Nx - 1)  # Grid spacing in the x direction, m
dy = dx
dt = dx**2 / (4 * PDMS_thermal_diffusivity_m2ps)  # Time step, s
x = np.linspace(0, square_length_m, Nx)
y = x
X, Y = np.meshgrid(x, y)

## Heat source parameters ##
beam_radius_m = 0.1
loading = 1e-6  # mass fraction of CB in PDMS, g/g
abs_coeff = 1e+6  # absorption coefficient of CB in PDMS, m^-1

radius = (X - square_length_m / 2)**2 + (Y - square_length_m / 2)**2
beam_mask = radius <= beam_radius_m
q = np.zeros((Nx, Ny))
q[beam_mask] = Q / (np.pi * beam_radius_m**2 * square_length_m)  # power density distributed across the entire volume
q *= np.exp(-1 * abs_coeff * loading * Y)  # depth-dependent exponential decay

# Show what fraction of power extends beyond the cube
transmittance = (Q - np.sum(q * dx**2)) / Q
print(f'Transmittance through material: {transmittance:.2%}')
print(f'max q (W): {np.max(q) * dx**2:.2e}')

# Plotting q's distribution
q_center = q[Nx//2, Ny//2, :]
plt.figure(figsize=(8, 6))
plt.plot(z, q_center)
plt.xlabel('Depth (m)')
plt.ylabel('Power Density (W/m^3)')
plt.title('Power Density Distribution along Beam Central Axis')
plt.grid(True)
plt.show()

# q = 10

def Runge_Kutta(T):
    """computes T at the the next time step with a 4th degree Runge-Kutta"""
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
        
    return T + (k1 + 2*k2 + 2*k3 + k4) / 6

def Boundary_Conditions(T_new):
    '''applies boundary conditions and surface convection to T_new'''
    ## boundary conditions ##
    T_new[0, :, :] = T_edge
    T_new[-1, :, :] = T_edge
    T_new[:, 0, :] = T_edge
    T_new[:, -1, :] = T_edge
    T_new[:, :, 0] = T_edge
    T_new[:, :, -1] = T_new[:, :, -1] - h_conv * (T_new[:, :, -1] - T_air) * dt / (dz) # convection at top surface
    return T_new

def Compute_T(output_times):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    # initialize temperature distribution and output structures
    T = np.full((Nx, Ny, Nz), T_edge)
    output_indices = [int(t / dt) for t in output_times] # times enumerated as indices based on time step
    output_temperatures = []

    # progress bar
    total_iterations = max(output_indices) + 1
    one_percent_iterations = int(total_iterations * 0.01)
    pbar = tqdm(total=total_iterations, desc="Computing", postfix="0% done")

    # loop across each necessary time step
    for n in range(max(output_indices) + 1):
        T_new = Runge_Kutta(T)
        T = Boundary_Conditions(T_new)
        
        if n in output_indices:
            output_temperatures.append(T.copy())

        if np.any(T) > 26: # troubleshooting problem where T never changes but the plots show it does
            print(T)
           
        # if n % one_percent_iterations == 0:
        #     progress_percentage = n * 100 / total_iterations
        #     pbar.set_postfix_str(f"{progress_percentage:.2f}% done")
        #     pbar.update(one_percent_iterations)

    pbar.close()
    return output_temperatures

def Plot_T_Slices(X, Y, Z, output_temperatures, output_times):
    '''plots slices of T at each output time'''
    fig = plt.figure(figsize=(12, 8))
    custom_cmap = plt.cm.get_cmap('hot', 20) # also, 'coolwarm'
    for i, T_out in enumerate(output_temperatures):
        ax = fig.add_subplot(1, len(output_times), i+1, projection='3d')
        c = ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, Nz//2], 
                        levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, Nz-1], 
                        levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, 0], 
                        levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_zlabel('z (m)')
        ax.set_xlim([0, square_length_m])
        ax.set_ylim([0, square_length_m])
        ax.set_zlim([0, square_length_m])
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=np.arange(20, 116, 10))
    plt.tight_layout()
    plt.show()

def save_output_temperatures(output_temperatures, filename):
    """
    Save the output temperatures to a file.
    
    Parameters:
    - output_temperatures: List of 3D numpy arrays containing temperature data.
    - filename: Name of the file to save the data to.
    """
    np.savez(filename, *output_temperatures)

output_times = [0, 1, 3, 5]
output_temperatures_RK = Compute_T(output_times)

## plot and/or save to file as npz ##
Plot_T_Slices(X, Y, Z, output_temperatures_RK, output_times)
save_output_temperatures(output_temperatures_RK, "output_temperatures.npz")